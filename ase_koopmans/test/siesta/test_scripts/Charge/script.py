# In this script the Virtual Crystal approximation is used to model
# a stronger affinity for positive charge on the H atoms.
# This could model interaction with other molecules not explicitly
# handled.
import numpy as np
from ase_koopmans.calculators.siesta import Siesta
from ase_koopmans.calculators.siesta.parameters import Species
from ase_koopmans.optimize import QuasiNewton
from ase_koopmans import Atoms

atoms = Atoms('CH4', np.array([
    [0.000000, 0.000000, 0.000000],
    [0.682793, 0.682793, 0.682793],
    [-0.682793, -0.682793, 0.682790],
    [-0.682793, 0.682793, -0.682793],
    [0.682793, -0.682793, -0.682793]]))

siesta = Siesta(
    species=[
        Species(symbol='H', excess_charge=0.1)])

atoms.calc = siesta
dyn = QuasiNewton(atoms, trajectory='h.traj')
dyn.run(fmax=0.02)
