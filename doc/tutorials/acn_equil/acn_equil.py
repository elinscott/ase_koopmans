import numpy as np

import ase_koopmans.units as units
from ase_koopmans import Atoms
from ase_koopmans.calculators.acn import ACN, m_me, r_cn, r_mec
from ase_koopmans.constraints import FixLinearTriatomic
from ase_koopmans.io import Trajectory
from ase_koopmans.md import Langevin

pos = [[0, 0, -r_mec],
       [0, 0, 0],
       [0, 0, r_cn]]
atoms = Atoms('CCN', positions=pos)
atoms.rotate(30, 'x')

# First C of each molecule needs to have the mass of a methyl group
masses = atoms.get_masses()
masses[::3] = m_me
atoms.set_masses(masses)

# Determine side length of a box with the density of acetonitrile at 298 K
d = 0.776 / 1e24  # Density in g/Ang3 (https://pubs.acs.org/doi/10.1021/je00001a006)
L = ((masses.sum() / units.mol) / d)**(1 / 3.)
# Set up box of 27 acetonitrile molecules
atoms.set_cell((L, L, L))
atoms.center()
atoms = atoms.repeat((3, 3, 3))
atoms.set_pbc(True)

# Set constraints for rigid triatomic molecules
nm = 27
atoms.constraints = FixLinearTriatomic(
    triples=[(3 * i, 3 * i + 1, 3 * i + 2)
             for i in range(nm)])

tag = 'acn_27mol_300K'
atoms.calc = ACN(rc=np.min(np.diag(atoms.cell)) / 2)

# Create Langevin object
md = Langevin(atoms, 1 * units.fs,
              temperature=300 * units.kB,
              friction=0.01,
              logfile=tag + '.log')

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(5000)

# Repeat box and equilibrate further
atoms.set_constraint()
atoms = atoms.repeat((2, 2, 2))
nm = 216
atoms.constraints = FixLinearTriatomic(
    triples=[(3 * i, 3 * i + 1, 3 * i + 2)
             for i in range(nm)])

tag = 'acn_216mol_300K'
atoms.calc = ACN(rc=np.min(np.diag(atoms.cell)) / 2)

# Create Langevin object
md = Langevin(atoms, 2 * units.fs,
              temperature=300 * units.kB,
              friction=0.01,
              logfile=tag + '.log')

traj = Trajectory(tag + '.traj', 'w', atoms)
md.attach(traj.write, interval=1)
md.run(3000)
