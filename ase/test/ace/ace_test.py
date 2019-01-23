import sys
from ase import Atoms
from ase.calculators.ace_cal import ACE

label = sys.argv[1]    
#mol= read('H2.xyz')
mol = Atoms('H2',[(0, 0, 0),(0, 0, 0.7)])
ace = ACE(label=label,command = 'mpirun -np 4 /home/khs/hs_file/programs/ACE-Molecule/ace PREFIX.inp > PREFIX.log')
mol.set_calculator(ace)
print ( mol.get_potential_energy())
