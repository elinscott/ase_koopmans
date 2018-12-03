from ase import Atoms
from ase.constraints import FixLinearTriatomic
from ase.calculators.acn import (ACN, m_me, m_c,
                                 m_n, r_mec, r_cn)
from ase.md import Langevin
import ase.units as units
from ase.io import Trajectory

import numpy as np


# Set up acetonitrile box at 25 deg C density
pos = [[0, 0, -r_mec],
       [0, 0, 0],
       [0, 0, r_cn]]
atoms = Atoms('CCN', positions=pos)
atoms.rotate(30, 'x')

# First C of each molecule needs to have the mass of a methyl group
masses = atoms.get_masses()
masses[::3] = m_me
atoms.set_masses(masses)

vol = ((masses.sum() / 6.022140857e23) / (0.776 / 1e24))**(1 / 3.)
atoms.set_cell((vol, vol, vol))
atoms.center()

atoms = atoms.repeat((3, 3, 3))
atoms.set_pbc(True)

# Set constraints for rigid triatomic molecules   
nm = 27
atoms.constraints = FixLinearTriatomic(
                    pairs=[(3 * i, 3 * i + 2)
                           for i in range(nm)],
                    centers=[j * 3 + 1 for j in range(nm)],
                    distances=[r_mec,r_cn],
                    masses=[m_me,m_c,m_n])

tag = 'acn_27mol_300K'
atoms.calc = ACN(rc=np.min(np.diag(atoms.cell))/2)

# Create Langevin object 
md = Langevin(atoms, 1 * units.fs, 
              temperature=300 * units.kB,
              selectlinear=[],
              friction=0.01, 
              logfile=tag + '.log')

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(5000)

# Repeat box and equilibrate further
tag = 'acn_216mol_300K'
atoms.set_constraint()
atoms = atoms.repeat((2, 2, 2))
nm = 216
atoms.constraints = FixLinearTriatomic(
                    pairs=[(3 * i, 3 * i + 2)
                           for i in range(nm)],
                    centers=[j * 3 + 1 for j in range(nm)],
                    distances=[r_mec,r_cn],
                    masses=[m_me,m_c,m_n])
atoms.calc = ACN(rc=np.min(np.diag(atoms.cell))/2)

# Create Langevin object
md = Langevin(atoms, 2 * units.fs,
              temperature=300 * units.kB,
              selectlinear=[],
              friction=0.01,
              logfile=tag + '.log')

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(3000)

