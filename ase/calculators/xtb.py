import os 
import numpy as np
from warnings import warn
from ase.calculators.calculator import FileIOCalculator
from ase.io.xyz import simple_write_xyz
from ase.units import Hartree, Bohr

class XTB(FileIOCalculator):
    implemented_properties = ['energy', 'forces']

    name = 'xTB'
    command = 'xtb PREFIX.xyz > PREFIX.out'

    if not 'XTBHOME' in os.environ:  # program installed?
        warn('WARNING: $XTBHOME variable not found.\
              xTB installed?')

    def __init__(self, label='ase_xtb', atoms=None, charge=None,
                 unpaired_electrons=None, gbsa=None):
        self.label = label
        self.atoms = atoms
        self.charge = charge
        self.uhf = unpaired_electrons
        self.gbsa = gbsa 

        add_pfx = ''
        if charge:
            add_pfx += ' --chrg {} '.format(charge)
        if unpaired_electrons:
            add_pfx += ' --uhf {} '.format(charge)
        if gbsa:
            add_pfx += ' --gbsa {} '.format(gbsa)

        self.command = command.split(' ')[0] +\
                       add_pfx +\
                       command.split('xtb ')[-1]

        FileIOCalculator.__init__(self, restart=None, 
                                  ignore_bad_restart_file=False,
                                  label=label,
                                  atoms=atoms, 
                                  command=self.command)


    def write_input(self, atoms, properties=None, system_changes=None):
        if 'forces' in properties:
            new_cmd = self.command.split(' ')
            new_cmd = [new_cmd[0]] + ['--grad'] + new_cmd[1:]
            new_cmd = ' '.join(x for x in new_cmd)
            self.command = new_cmd

        fobj = open(self.label + '.xyz', 'w')
        simple_write_xyz(fobj, [atoms],
                         'input xyz for xTB calc written by ASE')


    def read_results(self):
        self.read_energy()
        if '--grad' in self.command:
            self.read_forces()

    def read_energy(self):
        ''' xTB simply outputs to file: 'energy', so read from 
            that for now. We could parse the .out later for more stuff. 
            'energy' and 'forces' are appended to. Read last entry '''

        with open('energy', 'r') as f:
            lines = f.readlines()

        # always only 1 footerline
        energy = float(lines[-2].split()[-1])

        self.results['energy'] = energy * Hartree

    def read_forces(self):
        ''' file has 2 headerlines, then atomic positions, then forces.'''
        with open('gradient', 'r') as f:
            lines = f.readlines()

        # strip header, go to last cycle
        lines = lines[2:-1]
        cyc_idx = max([loc for loc, line in enumerate(lines) if 'cycle' in line])
        lines = lines[cyc_idx + 1:]
        idx = len(atoms)
        lines = lines[idx:]
        forces = np.zeros((idx, 3))
        for l, line in enumerate(lines):
            forces[l] = np.array([float(x.replace('D', 'e')) 
                                 for x in line.split()])

        forces *= Hartree / Bohr
        self.results['forces'] = forces


        
