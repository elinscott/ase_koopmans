"""This module defines an ASE interface to ORCA 4
by Ragnar Bjornsson
Based on NWchem interface but simplified.
Only supports energies and gradients (no dipole moments, orbital energies etc.) for now.
For more ORCA-keyword flexibility, method/xc/basis etc. keywords are not used.
Instead two keywords, orcasimpleinput and orcablock are used to define
the ORCA simple-inputline and the ORCA-block input.
This allows for more flexible use of any ORCA method or keyword available in ORCA
instead of hardcoding stuff.

Point Charge IO functionality added by A. Dohn.
"""
import os
import numpy as np

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.orca import write_orca
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError


class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []


class ORCA(FileIOCalculator):
    implemented_properties = ['energy', 'forces']

    if 'ORCA_COMMAND' in os.environ:
        command = os.environ['ORCA_COMMAND'] + ' PREFIX.inp > PREFIX.out'
    else:
        command = 'orca PREFIX.inp > PREFIX.out'

    default_parameters = dict(
        charge=0, mult=1,
        task='gradient',
        orcasimpleinput='PBE def2-SVP',
        orcablocks='%scf maxiter 200 end',
        )  

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='orca', atoms=None, **kwargs):
        """Construct ORCA-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()


    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters

        if self.pcpot:  # also write point charge file and add things to input
            pcstring = '% pointcharges \"' +\
                        self.label +\
                        '.pc\"\n\n% method \nDoEQ true \nend\n\n' 
            p['orcablocks'] += pcstring
            self.pcpot.write_mmcharges(self.label)

        p.write(self.label + '.ase')
        f = open(self.label + '.inp', 'w')
        f.write("! tightscf engrad %s \n" % p.orcasimpleinput);
        f.write("%s \n" % p.orcablocks)

        write_orca(f, atoms, p.charge, p.mult)
        f.close()

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        f = open(self.label + '.inp')
        for line in f:
            if line.startswith('geometry'):
                break
        symbols = []
        positions = []
        for line in f:
            if line.startswith('end'):
                break
            words = line.split()
            symbols.append(words[0])
            positions.append([float(word) for word in words[1:]])

        self.parameters = Parameters.read(self.label + '.ase')
        self.read_results()

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()

    def read_energy(self):
        """Read Energy from ORCA output file."""
        text = open(self.label + '.out', 'r').read()
        lines = iter(text.split('\n'))
        # Energy:
        estring = 'FINAL SINGLE POINT ENERGY'
        for line in lines:
            if estring in line:
                energy = float(line.split()[-1])
                break
        self.results['energy'] = energy * Hartree

    def read_forces(self):
        """Read Forces from ORCA output file."""
        file = open('orca.engrad', 'r')
        lines = file.readlines()
        file.close()
        getgrad="no"
        for i, line in enumerate(lines):
            if line.find('# The current gradient') >= 0:
                getgrad="yes";
                gradients = []
                tempgrad=[]
                continue
            if getgrad=="yes" and "#" not in line:
                grad=line.split()[-1]
                tempgrad.append(float(grad))
                if len(tempgrad)==3:
                    gradients.append(tempgrad)
                    tempgrad=[]
            if '# The at' in line:
                getgrad="no"
        self.results['forces'] = -np.array(gradients) * Hartree / Bohr

    def embed(self, mmcharges=None, **parameters):
        """Embed atoms in point-charges (mmcharges)
        """
        self.pcpot = PointChargePotential(mmcharges)
        return self.pcpot
        

class PointChargePotential:
    def __init__(self, mmcharges, positions=None, directory=None):
        """ Point Charge Potential Interface to ORCA """
        if positions is not None:
            self.set_positions(positions)
        if directory is None:
            directory = os.getcwd()

        self.directory = directory + os.sep
        self.mmcharges = mmcharges

    def set_positions(self, positions):
        self.positions = positions

    def set_charges(self, mmcharges):
        self.q_p = mmcharges

    def write_mmcharges(self, filename='orca_mm'):
        pc_file = open(os.path.join(self.directory, 
                                    filename + '.pc'), 'w')

        pc_file.write('{0:d}\n'.format(len(self.mmcharges)))
        for [pos, pc] in zip(self.positions, self.mmcharges):
            [x, y, z] = pos
            pc_file.write('{0:12.6f} {1:12.6f} {2:12.6f} {3:12.6f}\n'
                          .format(pc, x, y, z))








