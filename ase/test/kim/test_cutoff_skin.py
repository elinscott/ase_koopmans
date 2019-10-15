"""
To test that the calculator handles skin and cutoffs correctly.  Skin
should be added to both the influence distance and cutoff.  The cutoff
for the model ex_model_Ar_P_Morse_07C is 8.15 Angstroms and the default
skin distance when using the kimpy neighbor list library (which is the
default) is 0.2 times the cutoff distance (.63 Angstroms for this
model).
"""
import numpy as np
from ase.calculators.kim import KIM
from ase import Atoms


# Create calculator
modelname = "ex_model_Ar_P_Morse_07C"
calc = KIM(modelname)

# Create an FCC argon crystal
atoms = Atoms("Ar2", positions=[[0, 0, 0], [0, 0, 8.2]])

# Attach calculator to the atoms
atoms.set_calculator(calc)

# Get energy and forces
e_ideal = atoms.get_potential_energy()

# Perturb the x coordinate of the first atom by less than the skin distance
atoms.positions[0, 2] += 0.1

# Get new energy
e_perturb = atoms.get_potential_energy()

assert not np.isclose(e_ideal, e_perturb)
