"""
To test that the calculator handles skin and cutoffs correctly.
Skin should be added to both the influence distance and cutoff.
The cutoff for the model ex_model_Ar_P_Morse_07C is 8.15
"""
import numpy as np
from ase.calculators.kim import KIM
from ase import Atoms


# create calculator
modelname = "ex_model_Ar_P_Morse_07C"
calc = KIM(modelname)

# create an FCC crystal
argon = Atoms("Ar2", positions=[[0, 0, 0], [0, 0, 8.2]])

# attach calculator to the atoms
argon.set_calculator(calc)

# get energy and forces
e0 = argon.get_potential_energy()

# perturb the x coord of the first atom
argon.positions[0, 2] += 0.1

# get energy and forces
e1 = argon.get_potential_energy()

assert not np.isclose(e0, e1)
