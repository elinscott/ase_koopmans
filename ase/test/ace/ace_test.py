import sys
import os
from ase.atoms import Atoms
from ase.io import read
from ase.calculators.ace_cal import ACE

label = sys.argv[1]    
mol= read('H2.xyz')
ace = ACE(label=label)
mol.set_calculator(ace)
print ( mol.get_potential_energy())
