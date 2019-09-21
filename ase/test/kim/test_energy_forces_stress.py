"""
To test that the calculator can produce correct energy and forces.
"""
import numpy as np
from ase.calculators.kim import KIM
from ase.lattice.cubic import FaceCenteredCubic


energy_ref = 19.7196709065
forces_ref = [[-0.33209865, -13.98929271, -13.98929271],
              [0.18090261, 13.9896848, -13.98691618],
              [0.18090261, -13.98691618, 13.9896848],
              [-0.02970657, 13.98652409, 13.98652409]]
stress_ref = [-2.21148294e+00, -1.55423383e+00, -1.55423383e+00,
              2.17827079e-05, -8.39978013e-03, -8.39978013e-03]


def test_main():
    # create calculator
    modelname = 'ex_model_Ar_P_Morse_07C'
    calc = KIM(modelname)

    # create an FCC crystal
    argon = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              size=(1, 1, 1), symbol='Ar', pbc=(1, 0, 0),
                              latticeconstant=3.0)

    # perturb the x coord of the first atom
    argon.positions[0, 0] += 0.01

    # attach calculator to the atoms
    argon.set_calculator(calc)

    # get energy and forces
    energy = argon.get_potential_energy()
    forces = argon.get_forces()
    stress = argon.get_stress()

    tol = 1e-6
    assert np.isclose(energy, energy_ref, tol)
    assert np.allclose(forces, forces_ref, tol)
    assert np.allclose(stress, stress_ref, tol)

    # This has been known to segfault
    argon.set_pbc(True)
    argon.get_potential_energy()


if __name__ == '__main__':
    test_main()
