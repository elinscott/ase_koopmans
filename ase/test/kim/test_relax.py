"""Test that KIM works with a relaxation"""

import numpy as np
from ase.cluster import Icosahedron
from ase.calculators.kim import KIM
from ase.optimize import BFGS

energy_ref = -0.5420939378624228
positions_ref = np.array([[1.80450287, 1.80450287, 1.80450287],
                          [4.79991686, 1.80450287, -0.0467648],
                          [4.79991686, 1.80450287, 3.65577053],
                          [-1.19091113, 1.80450287, -0.0467648],
                          [-1.19091113, 1.80450287, 3.65577053],
                          [-0.0467648, 4.79991686, 1.80450287],
                          [3.65577053, 4.79991686, 1.80450287],
                          [-0.0467648, -1.19091113, 1.80450287],
                          [3.65577053, -1.19091113, 1.80450287],
                          [1.80450287, -0.0467648, 4.79991686],
                          [1.80450287, 3.65577053, 4.79991686],
                          [1.80450287, -0.0467648, -1.19091113],
                          [1.80450287, 3.65577053, -1.19091113]])


def test_relax():
    # Create structure
    atoms = Icosahedron('Ar', latticeconstant=3., noshells=2)

    # create calculator
    modelname = 'ex_model_Ar_P_Morse_07C'
    calc = KIM(modelname)

    # attach calculator to the atoms
    atoms.set_calculator(calc)

    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=0.05)

    assert np.isclose(atoms.get_potential_energy(), energy_ref)
    assert np.allclose(atoms.get_positions(), positions_ref)


if __name__ == '__main__':
    test_relax()
