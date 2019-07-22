"""
authors: Ben Comer (Georgia Tech), Xiangyun (Ray) Lei (Georgia Tech)

"""
from ase.calculators.calculator import Calculator, all_properties, all_changes
from ase.calculators.calculator import InputError, CalculationFailed, SCFError, ReadError 
import numpy as np
from ase.units import Bohr, Hartree
import warnings
import psi4
import os
import pickle
import codecs


class Psi4(Calculator):
    """
    An ase calculator for the popular open source Q-chem code
    psi4. This is really rudimentary

    you can always use the in-built psi4 module through:
    calc.psi4

    xc is the generic input for whatever method you wish to use, thus
    and quantum chemistry method implemented in psi4 can be input 
    (i.e. ccsd(t))
    """
    implemented_properties = ['energy', 'forces']
    
    default_parameters = {
                  "basis": "aug-cc-pvtz",
                  "num_threads": None,
                  "xc": "hf",
                  "memory": None,
                  'charge': None,
                  'multiplicity': None,
                  'reference': None,
                  'symmetry':'c1',
                  'PSI_SCRATCH' : '.',}
    def __init__(self, restart=None, ignore_bad_restart=False,
                 label='psi4-calc', atoms=None, command=None,
                 **kwargs):
        Calculator.__init__(self, restart=restart, ignore_bad_restart=ignore_bad_restart,
                            label=label, atoms=atoms, command=command,
                            **kwargs)
        self.psi4 = psi4
        # perform initial setup of psi4 python API
        self.set_psi4(atoms = atoms)

    def set_psi4(self, atoms = None):
        """
        This function sets the imported psi4 module to the settings the user 
        defines
        """
        if 'PSI_SCRATCH' in os.environ:
            pass
        else:
            os.environ['PSI_SCRATCH'] = self.parameters['PSI_SCRATCH']

        # Input spin settings
        if self.parameters['reference'] is not None:
            self.psi4.set_options({'reference':
                                   self.parameters['reference']})
        # Memory
        if self.parameters['memory'] is not None:
            self.psi4.set_memory(self.parameters['memory'])

        # Threads
        if self.parameters['num_threads'] == 'max':
            import multiprocessing
            self.psi4.set_num_threads(multiprocessing.cpu_count())
        elif  type(self.parameters['num_threads']) == int:
            self.psi4.set_num_threads(self.parameters['num_threads'])

        if self.parameters['xc'].lower() == 'lda':
            warnings.warn('Psi4 does not have LDA implemented, SVWN will be used instead')
            self.parameters['xc'] = 'svwn'

        # deal with some ASE specific inputs
        if 'kpts' in self.parameters:
            warnings.warn('psi4 is a non-periodic code, and thus does not '
                          'require k-points. This arguement will be ignored')

        if self.parameters['xc'].lower() == 'lda':
            warnings.warn('Psi4 does not have LDA implemented, SVWN will be used instead')
            self.parameters['xc'] = 'svwn'

        if 'nbands' in self.parameters:
            warnings.warn('psi4 does is a quantum chemistry program, and thus does '
                          'not operate on the basis of bands, please select a basis'
                          ' set instead. This input is ignored.')
        if 'smearing' in self.parameters:
            warnings.warn('Finite temperature DFT is not implemented in psi4 currently,'
                          ' thus a smearing argument cannot be utilized. This argument'
                          ' is ignored')

        if atoms is None:
            if self.atoms is None:
                return None
            else:
                atoms = self.atoms 
        if self.atoms is None:
            self.atoms = atoms
        # Generate the atomic input
        result = ''
        for atom in atoms:
            temp = '{}\t{:.15f}\t{:.15f}\t{:.15f}\n'.format(atom.symbol, \
            atom.position[0],atom.position[1], \
            atom.position[2])
            result += temp
        result += 'symmetry {}\n'.format(self.parameters['symmetry'])
        result += 'units angstrom\n'
        if self.parameters['charge'] is not None and \
                self.parameters['multiplicity'] is not None:
            result += '{} {}\n'.format(self.parameters['charge'],
                                     self.parameters['multiplicity'])
        elif self.parameters['charge'] is not None:
            warnings.warn('A charge was provided without a spin multiplicity.'
                          'A multiplicity of 1 is assumed')
            result += '{} 1\n'.format(self.parameters['charge'])

        elif self.parameters['multiplicity'] is not None:
            result += '0 {}\n'.format(self.parameters['multiplicity'])
        

        if not os.path.isdir(self.directory):
            os.mkdir(self.directory)
        self.psi4.core.set_output_file(self.label + '.dat',
                                       False)
        self.molecule = psi4.geometry(result)

    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def check_state(self, atoms):
        system_changes = Calculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        # There's no periodicity in psi4
        if 'cell' in system_changes:
            system_changes.remove('cell')
        return system_changes
    
    def read(self, label):
        """Read psi4 outputs made from this ASE calculator
        """
        filename = label + '.dat'
        if not os.path.isfile(filename):
            raise ReadError('Could not find the psi4 output file: ' + filename)

        f = open(filename, 'r')
        txt = f.read()
        f.close()

        if '!ASE Information\n' not in txt:
            raise Exception('the output file must be made by the ase psi4 '
                            'interface to be read by it')
        info = txt.split('!ASE Information\n')[1]
        saved_dict = pickle.loads(codecs.decode(info.encode(), "base64"))
        self.atoms = saved_dict['atoms']
        self.parameters = saved_dict['parameters']
        self.results = saved_dict['results']


    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes, symmetry = 'c1'):
        """Do the calculation.

        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.  Calculated properties should
        be inserted into results dictionary like shown in this dummy
        example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}

        """
        Calculator.calculate(self, atoms = atoms)
        if atoms == None:
            if self.atoms is None:
                raise InputError('An atoms object must be provided to perform a calculation')
            else:
                atoms = self.atoms
        elif self.atoms == None:
            self.atoms = atoms
        if atoms.get_initial_magnetic_moments().any() != 0:
            self.parameters['reference'] = 'uhf'
            self.parameters['multiplicity'] = None
        # this inputs all the settings into psi4
        self.set_psi4(atoms = atoms)

        # Set up the method
        method = self.parameters['xc']
        basis = self.parameters['basis']

        # Do the calculations
        for item in properties:
            if item == 'energy':
                energy = self.psi4.energy('{}/{}'.format(method,basis), 
                                      molecule = self.molecule,)
                # convert to eV
                self.results['energy'] = energy * Hartree
            if item == 'forces':
                grad, wf = self.psi4.driver.gradient('{}/{}'.format(method,basis),
                                                 return_wfn=True)
                # energy comes for free
                energy = wf.energy()
                self.results['energy'] = energy * Hartree
                # convert to eV/A
                # also note that the gradient is -1 * forces
                self.results['forces'] = -1 * np.array(grad) * Hartree / Bohr
        # dump the calculator info to the psi4 file
        save_atoms = self.atoms.copy()
        del save_atoms.calc
        save_dict = {'parameters': self.parameters,
                     'results': self.results,
                     'atoms': save_atoms}
        pickled = codecs.encode(pickle.dumps(save_dict), "base64").decode()
        self.psi4.core.print_out('!ASE Information\n')
        self.psi4.core.print_out(pickled)
