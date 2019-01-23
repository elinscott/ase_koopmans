from ase import Atoms
from ase.calculators.ace_cal import ACE

label = "test"
mol = Atoms('H2',[(0, 0, 0),(0, 0, 0.7)])
ace = ACE(label=label,command = '/home/khs/hs_file/programs/ACE-Molecule/ace PREFIX.inp > PREFIX.log')
mol.set_calculator(ace)
mol.get_potential_energy()
