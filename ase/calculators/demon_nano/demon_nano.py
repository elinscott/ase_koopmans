from __future__ import print_function
"""This module defines an ASE interface to deMon-nano.

http://website_deMon_nano

"""
import os
import os.path as op
import subprocess
import pickle
#import shutil

import numpy as np

from ase.units import Bohr, Hartree
import ase.data
from ase.calculators.calculator import FileIOCalculator, ReadError
from ase.calculators.calculator import Parameters, all_changes
from ase.calculators.calculator import equal
import ase.io

m_e_to_amu = 1822.88839

class Parameters_deMon_nano(Parameters):
    """Parameters class for the calculator.
    Documented in Base_deMon.__init__

    The options here are the most important ones that the user needs to be
    aware of. Further options accepted by deMon can be set in the dictionary
    input_arguments.

    """
    def __init__(
            self,
            label='rundir',
            atoms=None,
            command=None,
            restart=None,
            basis_path=None,
            ignore_bad_restart_file=False,
            deMon_restart_path='.',
            print_out='ASE',
            title='deMon input file',
            forces=False,
            basis={},
            input_arguments=None):
        kwargs = locals()
        kwargs.pop('self')
        Parameters.__init__(self, **kwargs)


class Demon_Nano(FileIOCalculator):
    """Calculator interface to the deMon-nano code. """

    implemented_properties = (
        'energy',
        'forces')

    def __init__(self, **kwargs):
        """ASE interface to the deMon-nano code.
        
        The deMon-nano code can be obtained from http://website_demonnano

        The DEMON_NANO_COMMAND environment variable must be set to run the executable, in bash it would be set along the lines of
        export DEMON_NANO_COMMAND="pathway-to-deMon-binary/deMon.username.x"

        Parameters:

        label : str 
            relative path to the run directory
        atoms : Atoms object
            the atoms object
        command  : str
            Command to run deMon. If not present the environment varable DEMON_NANO_COMMAND will be used
        restart  : str
            Relative path to ASE restart directory for parameters and atoms object and results
        basis_path  : str 
            Relative path to the directory containing BASIS, AUXIS, ECPS, MCPS and AUGMENT
        ignore_bad_restart_file : bool 
            Ignore broken or missing ASE restart files
            By default, it is an error if the restart
            file is missing or broken.
        deMon_restart_path  : str 
            Relative path to the deMon restart dir
        title : str 
            Title in the deMon input file.
        forces : bool
            If True a force calcilation is enforced
        print_out : str | list 
            Options for the printing in deMon
        basis : dict 
            Definition of basis sets.
        input_arguments : dict 
            Explicitly given input arguments. The key is the input keyword
            and the value is either a str, a list of str (will be written on the same line as the keyword),
            or a list of lists of str (first list is written on the first line, the others on following lines.)
        """
        
        parameters = Parameters_deMon_nano(**kwargs)
        
        # Setup the run command
        command = parameters['command']
        if command is None:
            command = os.environ.get('DEMON_NANO_COMMAND')

        if command is None:
            mess = 'The "DEMON_NANO_COMMAND" environment is not defined.'
            raise ValueError(mess)
        else:
            parameters['command'] = command
            
        # Call the base class.
        FileIOCalculator.__init__(
            self,
            **parameters)

    def __getitem__(self, key):
        """Convenience method to retrieve a parameter as
        calculator[key] rather than calculator.parameters[key]

            Parameters:
                key       : str, the name of the parameters to get.
        """
        return self.parameters[key]

    def set(self, **kwargs):
        """Set all parameters.

        Parameters:
            kwargs  : Dictionary containing the keywords for deMon
        """
        # Put in the default arguments.
        kwargs = self.default_parameters.__class__(**kwargs)

        if 'parameters' in kwargs:
            filename = kwargs.pop('parameters')
            parameters = Parameters.read(filename)
            parameters.update(kwargs)
            kwargs = parameters

        changed_parameters = {}

        for key, value in kwargs.items():
            oldvalue = self.parameters.get(key)
            if key not in self.parameters or not equal(value, oldvalue):
                changed_parameters[key] = value
                self.parameters[key] = value

        return changed_parameters

    def calculate(self,
                  atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        """Capture the RuntimeError from FileIOCalculator.calculate
        and add a little debug information from the deMon output.

        See base FileIocalculator for documentation.
        """

        if atoms is not None:
            self.atoms = atoms.copy()

        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please set $%s environment variable ' %
                               ('DEMON_NANO_COMMAND') +
                               'or supply the command keyword')
        command = self.command  
        olddir = os.getcwd()

        # basis path
        basis_path = self.parameters['basis_path']
        if basis_path is None:
            basis_path = os.environ.get('DEMON_BASIS_PATH')

        # uncomment this line for master vers
        #if basis_path is None:
        #    raise RuntimeError('Please set basis_path keyword,' +
        #                       ' or the DEMON_BASIS_PATH' +
        #                       ' environment variable')

        # go to directory and run calculation
        os.chdir(self.directory)
        errorcode = subprocess.call(command, shell=True)
            
        os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))

        try:
            self.read_results()
        except:
            with open(self.directory + '/deMon.out', 'r') as f:
                lines = f.readlines()
            debug_lines = 30
            print('##### %d last lines of the deMon.out' % debug_lines)
            for line in lines[-1*debug_lines:]:
                print(line.strip())
            print('##### end of deMon.out')
            raise RuntimeError

    def set_label(self, label):
        """Set label directory
        """
        self.label = label

    # in our case self.directory = self.label
        self.directory = self.label
        if self.directory == '':
            self.directory = os.curdir

    def write_input(self, atoms, properties=None, system_changes=None):
        """Write input (in)-file.
        See calculator.py for further details.
 
        Parameters:
             atoms        : The Atoms object to write.
             properties   : The properties which should be calculated.
             system_changes : List of properties changed since last run.
        
        """
        # Call base calculator.
        FileIOCalculator.write_input(
            self,
            atoms=atoms,
            properties=properties,
            system_changes=system_changes)
     
        if system_changes is None and properties is None:
            return
    
        filename = self.label + '/deMon.inp'

        add_print = ''

        # Start writing the file.
        with open(filename, 'w') as f:
            # write keyword argument keywords
            value = self.parameters['title']
            self._write_argument('TITLE', value, f)
            f.write('\n')

            # obtain forces through a single BOMD step
            # only if forces is in properties, or if keyword forces is True
            value = self.parameters['forces']
            if 'forces' in properties or value:
                self._write_argument('MDYNAMICS', 'ZERO', f)
                self._write_argument('MDSTEP', 'MAX=1', f)
                #default timestep is 0.25 fs id not enough - uncomment the line below
                #self._write_argument('TIMESTEP', '0.1', f)
                
                #add_print = add_print + ' ' + 'MD OPT'

            # print argument, here other options could change this
            value = self.parameters['print_out']
            assert(type(value) is str)
            value = value + add_print

            if not len(value) == 0:
                self._write_argument('PRINT', value, f)
                f.write('\n')

            # write general input arguments
            self._write_input_arguments(f)
            f.write('\n')

            # write geometry
            self._write_atomic_coordinates(f, atoms)

            # write pickle of Parameters
            #pickle.dump(self.parameters,
            #            open(self.label + '/deMon_parameters.pckl', 'wb'))

            # write xyz file for good measure.
            ase.io.write(self.label + '/deMon_atoms.xyz', self.atoms)
            
    def read(self, restart_path):
       """Read parameters from directory restart_path."""
    
       self.set_label(restart_path)

       if not op.exists(restart_path + '/deMon.inp'):
           raise ReadError('The restart_path file {0} does not exist'
                           .format(restart_path))
    
       if op.exists(restart_path + '/deMon_parameters.pckl'):
           parameters = pickle.load(open(restart_path +
                                     '/deMon_parameters.pckl', 'r'))
           self.parameters = parameters
     
       self.atoms = self.deMon_inp_to_atoms(restart_path + '/deMon.inp')
       
       self.read_results()
     
    def _write_input_arguments(self, f):
       """Write directly given input-arguments."""
       input_arguments = self.parameters['input_arguments']
    
       # Early return
       if input_arguments is None:
           return

       for key, value in input_arguments.items():
           self._write_argument(key, value, f)
    
    def _write_argument(self, key, value, f):
       """Write an argument to file.
       key :  a string coresponding to the input keyword
       value : the arguemnts, can be a string, a number or a list
       f :  and open file
       """
       if key == 'BASISPATH':    
       # Write a basis path to file.
       # Has to be in lowercase for deMon-nano to work
           line = value.lower()
           f.write(line)
           f.write('\n')
     
       elif not isinstance(value, (tuple, list)):
       # for only one argument, write on same line
           line = key.upper()
           line += ' ' + str(value).upper()
           f.write(line)
           f.write('\n')

       # for a list, write first argument on the first line,
       # then the rest on new lines
       else:
           line = key
           if not isinstance(value[0], (tuple, list)):
               for i in range(len(value)):
                   line += ' ' + str(value[i].upper())
               f.write(line)
               f.write('\n')
           else:
               for i in range(len(value)):
                   for j in range(len(value[i])):
                       line += ' ' + str(value[i][j]).upper()
                   f.write(line)
                   f.write('\n')
                   line = ''
                    
    def _write_atomic_coordinates(self, f, atoms):
        """Write atomic coordinates.
        Parameters:
        - f:     An open file object.
        - atoms: An atoms object.
        """
        #f.write('#\n')
        #f.write('# Atomic coordinates\n')
        #f.write('#\n')
        f.write('GEOMETRY CARTESIAN ANGSTROM\n')

        for i in range(len(atoms)):
            xyz = atoms.get_positions()[i]
            chem_symbol = atoms.get_chemical_symbols()[i]
            # the string below add numbers for atoms of the same type (H1,H2)
            #chem_symbol += str(i + 1)
             
            # the code below adds nuclear charge and atom masses to GEOMETRY 
            # if tag is set to 1 then we have a ghost atom,
            # set nuclear charge to 0
            #if(atoms.get_tags()[i] == 1):
            #    nuc_charge = str(0)
            #else:
            #    nuc_charge = str(atoms.get_atomic_numbers()[i])
            #mass = atoms.get_masses()[i]
               
            line = '{0:9s}'.format(chem_symbol).rjust(10) + ' '
            line += '{0:.5f}'.format(xyz[0]).rjust(10) + ' '
            line += '{0:.5f}'.format(xyz[1]).rjust(10) + ' '
            line += '{0:.5f}'.format(xyz[2]).rjust(10) + ' '
            #line += '{0:5s}'.format(nuc_charge).rjust(10) + ' '
            #line += '{0:.5f}'.format(mass).rjust(10) + ' '
           
            f.write(line)
            f.write('\n')

        f.write('\n')

# Analysis routines
    def read_results(self):
       """Read the results from output files."""
       self.read_energy()
       self.read_forces(self.atoms)
       #self.read_eigenvalues()
       #self.read_dipole()
       #self.read_xray()
    
    def read_energy(self):
       """Read energy from deMon.ase output file."""

       filename = self.label + '/deMon.ase'

       if op.isfile(filename):
           with open(filename, 'r') as f:
               lines = f.readlines()
          
       for i in range(len(lines)):
            if lines[i].startswith(' DFTB total energy [Hartree]'):
                self.results['energy'] = float(lines[i+1])
                break
       else:
           raise RuntimeError


    def read_forces(self, atoms):
        """Read forces from the deMon.ase file."""

        natoms = len(atoms)
        filename = self.label + '/deMon.ase'

        if op.isfile(filename):
            with open(filename, 'r') as f:
                lines = f.readlines()

                # find line where the forces start
                flag_found = False
                for i in range(len(lines)):
                    if lines[i].rfind('DFTB gradients at 0 time step in a.u.') > -1:
                        start = i + 1
                        flag_found = True
                        break

                if flag_found:
                    self.results['forces'] = np.zeros((natoms, 3), float)
                    for i in range(natoms):
                        line = [s for s in lines[i + start].strip().split(' ')
                                if len(s) > 0]
                        # forces with or without -1 ?! in deMon2k -1 is used
                        f = -1.*np.array([float(x) for x in line[1:4]])
                        # output forces in a.u.
                        #self.results['forces'][i, :] = f
                        # output forces with real dimension
                        self.results['forces'][i, :] = f * (Hartree / Bohr)

    def deMon_inp_to_atoms(self, filename):
       """Routine to read deMon.inp and convert it to an atoms object."""

       with open(filename, 'r') as f:
           lines = f.readlines()

       # find line where geometry starts
       for i in range(len(lines)):
           if lines[i].rfind('GEOMETRY') > -1:
               if lines[i].rfind('ANGSTROM'):
                   coord_units = 'Ang'
               elif lines.rfind('Bohr'):
                   coord_units = 'Bohr'
               ii = i
               break

       chemical_symbols = []
       xyz = []
       atomic_numbers = []
       masses = []

       for i in range(ii + 1, len(lines)):
           try:
               line = lines[i].split()

               if(len(line) > 0):
                   for symbol in ase.data.chemical_symbols:
                       found = None
                       if line[0].upper().rfind(symbol.upper()) > -1:
                           found = symbol
                           break
                   
                       if found is not None:
                           chemical_symbols.append(found)
                       else:
                           break

                       xyz.append([float(line[1]), float(line[2]), float(line[3])])
           
               if len(line) > 4:
                   atomic_numbers.append(int(line[4]))
           
               if len(line) > 5:
                   masses.append(float(line[5]))

           except:
               raise RuntimeError

       if coord_units == 'Bohr':
           xyz = xyz * Bohr

       natoms = len(chemical_symbols)

       # set atoms object
       atoms = ase.Atoms(symbols=chemical_symbols, positions=xyz)

       # if atomic numbers were read in, set them
       if(len(atomic_numbers) == natoms):
           atoms.set_atomic_numbers(atomic_numbers)
           
       # if masses were read in, set them
       if(len(masses) == natoms):
           atoms.set_masses(masses)
       
       return atoms
