# flake8: noqa
"""This module defines an ASE interface to ABINIT.

http://www.abinit.org/
"""

import os
from glob import glob
from os.path import join

import numpy as np

from ase.units import Bohr, Hartree, fs
from ase.data import chemical_symbols
from ase.io.abinit import read_abinit
from ase.calculators.calculator import FileIOCalculator, Parameters, kpts2mp, \
    ReadError


def write_files_file(fd, prefix, ppp_list):
    fd.write('%s\n' % (prefix + '.in'))  # input
    fd.write('%s\n' % (prefix + '.txt'))  # output
    fd.write('%s\n' % (prefix + 'i'))  # input
    fd.write('%s\n' % (prefix + 'o'))  # output

    # XXX:
    # scratch files
    #scratch = self.scratch
    #if scratch is None:
    #    scratch = dir
    #if not os.path.exists(scratch):
    #    os.makedirs(scratch)
    #fd.write('%s\n' % (os.path.join(scratch, prefix + '.abinit')))
    fd.write('%s\n' % (prefix + '.abinit'))
    # Provide the psp files
    for ppp in ppp_list:
        fd.write('%s\n' % (ppp)) # psp file path


class Abinit(FileIOCalculator):
    """Class for doing ABINIT calculations.

    The default parameters are very close to those that the ABINIT
    Fortran code would use.  These are the exceptions::

      calc = Abinit(label='abinit', xc='LDA', ecut=400, toldfe=1e-5)
    """

    implemented_properties = ['energy', 'forces', 'stress', 'magmom']
    command = 'abinit < PREFIX.files > PREFIX.log'

    default_parameters = dict(
        xc='LDA',
        smearing=None,
        kpts=None,
        charge=0.0,
        raw=None,
        pps='fhi')

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='abinit', atoms=None, scratch=None, **kwargs):
        """Construct ABINIT-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'abinit'.

        Examples
        ========
        Use default values:

        >>> h = Atoms('H', calculator=Abinit(ecut=200, toldfe=0.001))
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()

        """

        self.scratch = scratch

        self.species = []
        self.ppp_list = []

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=tuple(), system_changes=tuple(),
        raise_exception=True):
        """Write input parameters to files-file."""

        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        try:
            if ('numbers' in system_changes or
                'initial_magmoms' in system_changes):
                self.initialize(atoms, raise_exception=raise_exception)
        except Exception as e:
            print(e, '... but I continue to complete the abinit.in')
            pass

        with open(self.label + '.files', 'w') as fd:
            write_files_file(fd, self.prefix, self.ppp_list)

        # Abinit will write to label.txtA if label.txt already exists,
        # so we remove it if it's there:
        filename = self.label + '.txt'
        if os.path.isfile(filename):
            os.remove(filename)

        param = self.parameters
        param.write(self.label + '.ase')

        from ase.io.abinit import write_abinit_in
        with open(self.label + '.in', 'w') as fd:
            write_abinit_in(fd, atoms, param, self.species)

    def read(self, label):
        """Read results from ABINIT's text-output file."""
        FileIOCalculator.read(self, label)
        filename = self.label + '.txt'
        if not os.path.isfile(filename):
            raise ReadError('ABINIT output file '+filename+' is missing.')

        self.atoms = read_abinit(self.label + '.in')
        self.parameters = Parameters.read(self.label + '.ase')

        self.initialize(self.atoms)
        self.read_results()

    def read_results(self):
        filename = self.label + '.txt'
        text = open(filename).read().lower()

        for line in iter(text.split('\n')):
            if line.rfind('error') > -1 or line.rfind('was not enough scf cycles to converge') > -1:
                raise ReadError(line)
            if line.rfind('natom  ') > -1:
                natoms = int(line.split()[-1])

        lines = iter(text.split('\n'))
        # Stress:
        # Printed in the output in the following format [Hartree/Bohr^3]:
        # sigma(1 1)=  4.02063464E-04  sigma(3 2)=  0.00000000E+00
        # sigma(2 2)=  4.02063464E-04  sigma(3 1)=  0.00000000E+00
        # sigma(3 3)=  4.02063464E-04  sigma(2 1)=  0.00000000E+00
        for line in lines:
            if line.rfind(
                'cartesian components of stress tensor (hartree/bohr^3)') > -1:
                stress = np.empty(6)
                for i in range(3):
                    entries = next(lines).split()
                    stress[i] = float(entries[2])
                    stress[i + 3] = float(entries[5])
                self.results['stress'] = stress * Hartree / Bohr**3
                break
        else:
            raise RuntimeError

        # Energy [Hartree]:
        # Warning: Etotal could mean both electronic energy and free energy!
        etotal = None
        efree = None
        if 'PAW method is used'.lower() in text:  # read DC energy according to M. Torrent
            for line in iter(text.split('\n')):
                if line.rfind('>>>>> internal e=') > -1:
                    etotal = float(line.split('=')[-1])*Hartree  # second occurrence!
            for line in iter(text.split('\n')):
                if line.rfind('>>>> etotal (dc)=') > -1:
                    efree = float(line.split('=')[-1])*Hartree
        else:
            for line in iter(text.split('\n')):
                if line.rfind('>>>>> internal e=') > -1:
                    etotal = float(line.split('=')[-1])*Hartree  # first occurrence!
                    break
            for line in iter(text.split('\n')):
                if line.rfind('>>>>>>>>> etotal=') > -1:
                    efree = float(line.split('=')[-1])*Hartree
        if efree is None:
            raise RuntimeError('Total energy not found')
        if etotal is None:
            etotal = efree

        # Energy extrapolated to zero Kelvin:
        self.results['energy'] = (etotal + efree) / 2
        self.results['free_energy'] = efree

        # Forces:
        for line in lines:
            if line.rfind('cartesian forces (ev/angstrom) at end:') > -1:
                forces = []
                for i in range(natoms):
                    forces.append(np.array(
                            [float(f) for f in next(lines).split()[1:]]))
                self.results['forces'] = np.array(forces)
                break
        else:
            raise RuntimeError
        #
        self.width = self.read_electronic_temperature()
        self.nband = self.read_number_of_bands()
        self.niter = self.read_number_of_iterations()
        self.nelect = self.read_number_of_electrons()
        self.results['magmom'] = self.read_magnetic_moment()


    def initialize(self, atoms, raise_exception=True):

        self.species = list(set(atoms.get_atomic_numbers()))
        self.spinpol = atoms.get_initial_magnetic_moments().any()
        self.ppp_list = self.get_ppp_list(self.species, atoms, raise_exception)


    def get_ppp_list(self, species, atoms, raise_exception):

        ppp_list = []

        pppaths = os.environ.get('ABINIT_PP_PATH','.').split(':')
        xcname = 'GGA' if self.parameters.xc != 'LDA' else 'LDA'
        pps = self.parameters.pps
        for Z in species:
            number = abs(Z)
            symbol = chemical_symbols[number]

            names  = []
            for s in [ symbol, symbol.lower() ]:
                for xcn in [ xcname, xcname.lower() ]:
                    if pps in ['paw']:
                        hghtemplate = '%s-%s-%s.paw'  # E.g. "H-GGA-hard-uspp.paw"
                        names.append(hghtemplate % (s, xcn, '*'))
                        names.append('%s[.-_]*.paw'   % s)
                    elif pps in ['pawxml']:
                        hghtemplate = '%s.%s%s.xml'  # E.g. "H.GGA_PBE-JTH.xml"
                        names.append(hghtemplate % (s, xcn, '*'))
                        names.append('%s[.-_]*.xml'   % s)
                    elif pps in ['hgh.k']:
                        hghtemplate = '%s-q%s.hgh.k'  # E.g. "Co-q17.hgh.k"
                        names.append(hghtemplate % (s, '*'))
                        names.append('%s[.-_]*.hgh.k' % s)
                        names.append('%s[.-_]*.hgh' % s)
                    elif pps in ['tm']:
                        hghtemplate = '%d%s%s.pspnc'  # E.g. "44ru.pspnc"
                        names.append(hghtemplate % (number, s, '*'))
                        names.append('%s[.-_]*.pspnc' % s)
                    elif pps in ['hgh', 'hgh.sc']:
                        hghtemplate = '%d%s.%s.hgh'  # E.g. "42mo.6.hgh"
                        # There might be multiple files with different valence
                        # electron counts, so we must choose between
                        # the ordinary and the semicore versions for some elements.
                        #
                        # Therefore we first use glob to get all relevant files,
                        # then pick the correct one afterwards.
                        names.append(hghtemplate % (number, s, '*'))
                        names.append('%d%s%s.hgh' % (number, s, '*'))
                        names.append('%s[.-_]*.hgh' % s)
                    else: # default extension
                        names.append('%02d-%s.%s.%s' % (number, s, xcn, pps))
                        names.append('%02d[.-_]%s*.%s'   % (number, s, pps))
                        names.append('%02d%s*.%s'   % (number, s, pps))
                        names.append('%s[.-_]*.%s'   % (s, pps))

            found = False
            for name in names:        # search for file names possibilities
                for path in pppaths:  # in all available directories
                    filenames = glob(join(path, name))
                    if not filenames:
                        continue
                    if pps == 'paw':
                        # warning: see download.sh in
                        # abinit-pseudopotentials*tar.gz for additional
                        # information!
                        filenames[0] = max(filenames)  # Semicore or hard
                    elif pps == 'hgh':
                        filenames[0] = min(filenames)  # Lowest valence electron count
                    elif pps == 'hgh.k':
                        filenames[0] = max(filenames)  # Semicore - highest electron count
                    elif pps == 'tm':
                        filenames[0] = max(filenames)  # Semicore - highest electron count
                    elif pps == 'hgh.sc':
                        filenames[0] = max(filenames)  # Semicore - highest electron count

                    if filenames:
                        found = True
                        ppp_list.append(filenames[0])
                        break
                if found:
                    break

            if not found:
                ppp_list.append("Provide {}.{}.{}?".format(symbol, '*', pps))
                if raise_exception:
                    raise RuntimeError('Could not find {} pseudopotential {} for {}'.format(xcname.lower(), pps, symbol))

        return ppp_list


    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = None
        for line in open(self.label + '.txt'):
            if line.find(' At SCF step') != -1: # find the last iteration number
                niter = int(line.split()[3].rstrip(','))
        return niter

    def get_electronic_temperature(self):
        return self.width * Hartree

    def read_electronic_temperature(self):
        width = None
        # only in log file!
        for line in open(self.label + '.log'):  # find last one
            if line.find('tsmear') != -1:
                width = float(line.split()[1].strip())
        return width

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelect = None
        # only in log file!
        for line in open(self.label + '.log'):  # find last one
            if line.find('with nelect') != -1:
                nelect = float(line.split('=')[1].strip())
        return nelect

    def get_number_of_bands(self):
        return self.nband

    def read_number_of_bands(self):
        nband = None
        for line in open(self.label + '.txt'): # find last one
            if line.find('     nband') != -1: # nband, or nband1, nband*
                nband = int(line.split()[-1].strip())
        return nband

    def get_kpts_info(self, kpt=0, spin=0, mode='eigenvalues'):
        return self.read_kpts_info(kpt, spin, mode)

    def get_k_point_weights(self):
        return self.get_kpts_info(kpt=0, spin=0, mode='k_point_weights')

    def get_bz_k_points(self):
        raise NotImplementedError

    def get_ibz_k_points(self):
        return self.get_kpts_info(kpt=0, spin=0, mode='ibz_k_points')

    def get_spin_polarized(self):
        return self.spinpol

    def get_number_of_spins(self):
        return 1 + int(self.spinpol)

    def read_magnetic_moment(self):
        magmom = None
        if not self.get_spin_polarized():
            magmom = 0.0
        else: # only for spinpolarized system Magnetisation is printed
            for line in open(self.label + '.txt'):
                if line.find('Magnetisation') != -1: # last one
                    magmom = float(line.split('=')[-1].strip())
        return magmom

    def get_fermi_level(self):
        return self.read_fermi()

    def get_eigenvalues(self, kpt=0, spin=0):
        return self.get_kpts_info(kpt, spin, 'eigenvalues')

    def get_occupations(self, kpt=0, spin=0):
        return self.get_kpts_info(kpt, spin, 'occupations')

    def read_fermi(self):
        """Method that reads Fermi energy in Hartree from the output file
        and returns it in eV"""
        E_f=None
        filename = self.label + '.txt'
        text = open(filename).read().lower()
        assert 'error' not in text
        for line in iter(text.split('\n')):
            if line.rfind('fermi (or homo) energy (hartree) =') > -1:
                E_f = float(line.split('=')[1].strip().split()[0])
        return E_f*Hartree

    def read_kpts_info(self, kpt=0, spin=0, mode='eigenvalues'):
        """ Returns list of last eigenvalues, occupations, kpts weights, or
        kpts coordinates for given kpt and spin.
        Due to the way of reading output the spins are exchanged in spin-polarized case.  """
        # output may look like this (or without occupation entries); 8 entries per line:
        #
        #  Eigenvalues (hartree) for nkpt=  20  k points:
        # kpt#   1, nband=  3, wtk=  0.01563, kpt=  0.0625  0.0625  0.0625 (reduced coord)
        #  -0.09911   0.15393   0.15393
        #      occupation numbers for kpt#   1
        #   2.00000   0.00000   0.00000
        # kpt#   2, nband=  3, wtk=  0.04688, kpt=  0.1875  0.0625  0.0625 (reduced coord)
        # ...
        #
        assert mode in ['eigenvalues', 'occupations', 'ibz_k_points',
                        'k_point_weights'], mode
        if self.get_spin_polarized():
            spin = {0: 1, 1: 0}[spin]
        if spin == 0:
            spinname = ''
        else:
            spinname = 'SPIN UP'.lower()
        # number of lines of eigenvalues/occupations for a kpt
        nband = self.get_number_of_bands()
        n_entries_float = 8  # float entries per line
        n_entry_lines = max(1, int((nband - 0.1) / n_entries_float) + 1)

        filename = self.label + '.txt'
        text = open(filename).read().lower()
        assert 'error' not in text
        lines = text.split('\n')
        text_list = []
        # find the beginning line of last eigenvalues
        contains_eigenvalues = 0
        for n, line in enumerate(lines):
            if spin == 0:
                if line.rfind('eigenvalues (hartree) for nkpt') > -1:
                #if line.rfind('eigenvalues (   ev  ) for nkpt') > -1: #MDTMP
                    contains_eigenvalues = n
            else:
                if (line.rfind('eigenvalues (hartree) for nkpt') > -1 and
                    line.rfind(spinname) > -1): # find the last 'SPIN UP'
                        contains_eigenvalues = n
        # find the end line of eigenvalues starting from contains_eigenvalues
        text_list = [lines[contains_eigenvalues]]
        for line in lines[contains_eigenvalues + 1:]:
            text_list.append(line)
            # find a blank line or eigenvalues of second spin
            if (not line.strip() or
                line.rfind('eigenvalues (hartree) for nkpt') > -1):
                break
        # remove last (blank) line
        text_list = text_list[:-1]

        assert contains_eigenvalues, 'No eigenvalues found in the output'

        n_kpts = int(text_list[0].split('nkpt=')[1].strip().split()[0])

        # get rid of the "eigenvalues line"
        text_list = text_list[1:]

        # join text eigenvalues description with eigenvalues
        # or occupation numbers for kpt# with occupations
        contains_occupations = False
        for line in text_list:
            if line.rfind('occupation numbers') > -1:
                contains_occupations = True
                break
        if mode == 'occupations':
            assert contains_occupations, 'No occupations found in the output'

        if contains_occupations:
            range_kpts = 2*n_kpts
        else:
            range_kpts = n_kpts

        values_list = []
        offset = 0
        for kpt_entry in range(range_kpts):
            full_line = ''
            for entry_line in range(n_entry_lines+1):
                full_line = full_line+str(text_list[offset+entry_line])
            first_line = text_list[offset]
            if mode == 'occupations':
                if first_line.rfind('occupation numbers') > -1:
                    # extract numbers
                    full_line = [float(v) for v in full_line.split('#')[1].strip().split()[1:]]
                    values_list.append(full_line)
            elif mode in ['eigenvalues', 'ibz_k_points', 'k_point_weights']:
                if first_line.rfind('reduced coord') > -1:
                    # extract numbers
                    if mode == 'eigenvalues':
                        full_line = [Hartree*float(v) for v in full_line.split(')')[1].strip().split()[:]]
                        #full_line = [float(v) for v in full_line.split(')')[1].strip().split()[:]] #MDTMP
                    elif mode == 'ibz_k_points':
                        full_line = [float(v) for v in full_line.split('kpt=')[1].strip().split('(')[0].split()]
                    else:
                        full_line = float(full_line.split('wtk=')[1].strip().split(',')[0].split()[0])
                    values_list.append(full_line)
            offset = offset+n_entry_lines+1

        if mode in ['occupations', 'eigenvalues']:
            return np.array(values_list[kpt])
        else:
            return np.array(values_list)
