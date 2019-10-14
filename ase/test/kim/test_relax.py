"""Test that KIM works with a relaxation"""
import numpy as np
from ase.cluster import Icosahedron
from ase.calculators.kim import KIM
from ase.optimize import BFGS

energy_ref = -0.5420939378624228

# Create structure
atoms = Icosahedron("Ar", latticeconstant=3.0, noshells=2)

# create calculator
modelname = "ex_model_Ar_P_Morse_07C"
calc = KIM(modelname)

# attach calculator to the atoms
atoms.set_calculator(calc)

opt = BFGS(atoms, logfile=None)
opt.run(fmax=0.05)

assert np.isclose(atoms.get_potential_energy(), energy_ref)
