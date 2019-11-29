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
                 label='abinit', atoms=None, **kwargs):
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

        species = list(set(atoms.numbers))
        ppp = self.get_ppp_list(species, atoms, raise_exception)
        with open(self.label + '.files', 'w') as fd:
            write_files_file(fd, self.prefix, ppp)

        # Abinit will write to label.txtA if label.txt already exists,
        # so we remove it if it's there:
        filename = self.label + '.txt'
        if os.path.isfile(filename):
            os.remove(filename)

        param = self.parameters
        param.write(self.label + '.ase')

        from ase.io.abinit import write_abinit_in
        with open(self.label + '.in', 'w') as fd:
            write_abinit_in(fd, atoms, param, species)

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
        from ase.io.abinit import read_abinit_out, read_abinit_log
        self.results = {}
        with open(filename) as fd:
            results = read_abinit_out(fd)
            self.results.update(results)
        with open(self.label + '.log') as fd:
            results = read_abinit_log(fd)
            self.results.update(results)

    def initialize(self, atoms, raise_exception=True):
        self.spinpol = atoms.get_initial_magnetic_moments().any()

    def get_ppp_list(self, species, atoms, raise_exception):
        return get_ppp_list(species, atoms, raise_exception,
                            xc=self.parameters.xc,
                            pps=self.parameters.pps)

    def get_number_of_iterations(self):
        return self.niter

    def get_electronic_temperature(self):
        return self.width * Hartree

    def get_number_of_electrons(self):
        return self.nelect

    def get_number_of_bands(self):
        return self.results['nbands']

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

    def get_fermi_level(self):
        return self.read_fermi()

    def get_eigenvalues(self, kpt=0, spin=0):
        return self.get_kpts_info(kpt, spin, 'eigenvalues')

    def get_occupations(self, kpt=0, spin=0):
        return self.get_kpts_info(kpt, spin, 'occupations')


def get_ppp_list(species, atoms, raise_exception, xc, pps):
    ppp_list = []

    pppaths = os.environ.get('ABINIT_PP_PATH', '.').split(':')
    xcname = 'GGA' if xc != 'LDA' else 'LDA'
    for Z in species:
        number = abs(Z)
        symbol = chemical_symbols[number]

        names  = []
        for s in [symbol, symbol.lower()]:
            for xcn in [xcname, xcname.lower()]:
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
