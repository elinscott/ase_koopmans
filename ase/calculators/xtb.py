import os 
import numpy as np
from warnings import warn
from ase.calculators.calculator import FileIOCalculator
from ase.io.xyz import simple_write_xyz
from ase.units import Hartree, Bohr

class XTB(FileIOCalculator):
    ''' XTB Interface for ASE 
        WIPs: 
              1. Has to calculate forces every time. (easy fix)
              2. Should use default_parameters dict instead. 
              
                                Asmus O. Dohn July 2019''' 

                                    
    implemented_properties = ['energy', 'forces']

    name = 'xTB'
    command = 'xtb PREFIX.xyz > PREFIX.out'

    if not 'XTBHOME' in os.environ:  # program installed?
        warn('WARNING: $XTBHOME variable not found.\
              xTB installed?')

    def __init__(self, label='ase_xtb', atoms=None, charge=None,
                 unpaired_electrons=None, gbsa=None, restart=True, procs=1,
                 acc=None):

        self.label = label
        self.atoms = atoms
        self.charge = charge
        self.uhf = unpaired_electrons
        self.gbsa = gbsa
        self.restart = restart
        self.procs = procs
        self.acc = acc

        command = 'xtb PREFIX.xyz > PREFIX.out'
        add_pfx = ' --grad '  #XXX WIP: currently forces always calculated.
        if charge:
            add_pfx += ' --chrg {} '.format(charge)
        if unpaired_electrons:
            add_pfx += ' --uhf {} '.format(charge)
        if not restart:
            add_pfx += ' --norestart '
        if procs != 1:
            #add_pfx += '-P {0:d} '.format(procs)
            os.environ['OMP_NUM_THREADS'] = '{0:d}'.format(procs)
        if acc is not None:
            add_pfx += ' --acc {0:2.4f} '.format(acc) 

        self.command = command.split(' ')[0] +\
                       add_pfx +\
                       command.split('xtb ')[-1]

        FileIOCalculator.__init__(self, restart=None, 
                                  ignore_bad_restart_file=False,
                                  label=label,
                                  atoms=atoms, 
                                  command=self.command)


    def write_input(self, atoms, properties=None, system_changes=None):
        #if 'forces' in properties:
        #    new_cmd = self.command.split(' ')
        #    new_cmd = [new_cmd[0]] + ['--grad'] + new_cmd[1:]
        #    new_cmd = ' '.join(x for x in new_cmd)
        #    self.command = new_cmd

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

        # XXX WIP: OOPS! energy file ONLY shows up if --grad is entered.
        # need to read from normal file as well
        with open('energy', 'r') as f:
            lines = f.readlines()

        # always only 1 footerline
        energy = float(lines[-2].split()[-1])

        self.results['energy'] = energy * Hartree

    def read_forces(self):
        ''' file has 2 headerlines, then atomic positions, then forces.'''
        with open('gradient', 'r') as f:
            lines = f.readlines()
        
        idx = len(self.atoms)
        # strip header and footer, go to last cycle
        #lines = lines[2:-1]
        cyc_idx = max([loc for loc, line in enumerate(lines) if 'cycle' in line])
        lines = lines[cyc_idx + 1:]
        ### test on positions and elements
        xtb_out = lines[:idx]
        xtb_pos = np.zeros((len(self.atoms), 3))
        xtb_sym = []
        for l, pos in enumerate(xtb_out):
            x, y, z, sym = pos.split()
            xtb_pos[l] = [float(x), float(y), float(z)]
            xtb_sym.append(sym.capitalize())
        
        xtb_pos *= Bohr 
        ase_pos = self.atoms.get_positions()
        # assert(abs(xtb_pos - ase_pos) < 1e-6).all(), 'ERROR in force readout: Positions'
        # XXX WIP: these positions are 
        ase_sym = self.atoms.get_chemical_symbols()
        assert(xtb_sym == ase_sym), 'ERROR in force readout: Atom sequence'

        lines = lines[idx:-1]  # remove $end footer
        forces = np.zeros((idx, 3))
        for l, line in enumerate(lines):
            forces[l] = np.array([float(x.replace('D', 'e')) 
                                 for x in line.split()])

        forces *= Hartree / Bohr
        forces *= -1
        self.results['forces'] = forces


        
