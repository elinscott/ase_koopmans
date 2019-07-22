import numpy as np
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import SCFError
import ase.units


class QChem(FileIOCalculator):
    """
    QChem calculator
    """
    name = 'QChem'

    implemented_properties = ['energy', 'forces']
    command = 'qchem PREFIX.inp PREFIX.out'

    # Following the minimal requirements given in
    # http://www.q-chem.com/qchem-website/manual/qchem43_manual/sect-METHOD.html
    default_parameters = {'method': 'hf',
                          'basis': '6-31G*',
                          'jobtype': 'force',
                          'charge': 0}

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='qchem', scratch=None, np=1, nt=1, basisfile=None,
                 ecpfile=None, pbs=False, atoms=None, **kwargs):
        """
        All options are copied verbatim to the input file.

        scratch: str
            path used for the scratch file
        np: int
            number of processors
        nt: int
            number of threads
        pbs: boolean
            command line flag for pbs scheduler (see QChem manual)
        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        # Augment the command by various flags
        if pbs:
            self.command = 'qchem -pbs '
        else:
            self.command = 'qchem '
        if np != 1:
            self.command += '-np %d ' % np
        if nt != 1:
            self.command += '-nt %d ' % nt
        self.command += 'PREFIX.inp PREFIX.out'
        if scratch is not None:
            self.command += ' %s' % scratch

        self.basisfile = basisfile
        self.ecpfile = ecpfile

    def read(self, label):
        raise NotImplementedError

    def read_results(self):
        filename = self.label + '.out'

        with open(filename, 'r') as fileobj:
            lineiter = iter(fileobj)
            for line in lineiter:
                if 'SCF failed to converge' in line:
                    raise SCFError()
                elif 'ERROR: alpha_min' in line:
                    # Even though it is not technically a SCFError:
                    raise SCFError()
                elif ' Total energy in the final basis set =' in line:
                    convert = ase.units.Hartree
                    self.results['energy'] = float(line.split()[8]) * convert
                elif ' Gradient of SCF Energy' in line:
                    iforces = []
                    # Skip first line containing atom numbering
                    next(lineiter)
                    while True:
                        # Get next line and cut off the component numbering and
                        # remove trailing characters ('\n' and stuff)
                        line = next(lineiter)[5:].rstrip()
                        # Cut in chunks of 12 symbols and convert into strings
                        # This is prefered over string.split() as the fields
                        # may overlap when the gradient gets large
                        Fx = list(map(
                            float,
                            [line[i:i + 12] for i in range(0, len(line), 12)]))
                        # Repeat for Fy and Fz
                        line = next(lineiter)[5:].rstrip()
                        Fy = list(map(
                            float,
                            [line[i:i + 12] for i in range(0, len(line), 12)]))
                        line = next(lineiter)[5:].rstrip()
                        Fz = list(map(
                            float,
                            [line[i:i + 12] for i in range(0, len(line), 12)]))
                        iforces.extend(zip(Fx, Fy, Fz))

                        # After three force components we expect either a
                        # separator line, which we want to skip, or the end of
                        # the gradient matrix which is characterized by the
                        # line ' Max gradient component'.
                        # Maybe change stopping criterion to be independent of
                        # next line. Eg. if not lineiter.next().startswith(' ')
                        if ' Max gradient component' in next(lineiter):
                            # Minus to convert from gradient to force
                            self.results['forces'] = np.array(iforces) * (
                                -ase.units.Hartree / ase.units.Bohr)
                            break

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        filename = self.label + '.inp'

        with open(filename, 'w') as fileobj:
            fileobj.write('$comment\n   ASE generated input file\n$end\n\n')

            fileobj.write('$rem\n')
            for prm in self.parameters:
                if prm not in ['charge', 'multiplicity']:
                    fileobj.write('   %-25s   %s\n' % (prm,
                                                       self.parameters[prm]))
            # Not even a parameters as this is an absolute necessity
            fileobj.write('   %-25s   %s\n' % ('sym_ignore', 'true'))
            fileobj.write('$end\n\n')

            fileobj.write('$molecule\n')
            # Following the example set by the gaussian calculator
            if ('multiplicity' not in self.parameters):
                tot_magmom = atoms.get_initial_magnetic_moments().sum()
                mult = tot_magmom + 1
            else:
                mult = self.parameters['multiplicity']
            # Default charge of 0 is defined in default_parameters
            fileobj.write('   %d %d\n' % (self.parameters['charge'], mult))
            for a in atoms:
                fileobj.write('   %s  %f  %f  %f\n' % (a.symbol,
                                                       a.x, a.y, a.z))
            fileobj.write('$end\n\n')

            if self.basisfile is not None:
                with open(self.basisfile, 'r') as f_in:
                    basis = f_in.readlines()
                fileobj.write('$basis\n')
                fileobj.writelines(basis)
                fileobj.write('$end\n\n')

            if self.ecpfile is not None:
                with open(self.ecpfile, 'r') as f_in:
                    ecp = f_in.readlines()
                fileobj.write('$ecp\n')
                fileobj.writelines(ecp)
                fileobj.write('$end\n\n')
