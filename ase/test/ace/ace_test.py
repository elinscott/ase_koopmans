import sys
import os
from ase.atoms import Atoms
from ase.io import read
from ase.calculators.ace_cal import ACE

Label = str(sys.argv[1].split('.inp')[0])    
mol= read('H2.xyz')
ace = ACE(label=Label)
mol.set_calculator(ace)
print ('Single Point Energy : ',  ace.get_property('energy', mol))
