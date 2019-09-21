import numpy as np
import ase.build
from ase.calculators.kim import KIM

energy_ref = -49.105522702876016
forces_ref = np.array([[5.27355937e-15, -1.97733646e-01, -4.60574165e-21],
                       [8.71264866e-16, -9.07596236e-03, 1.85559638e-01],
                       [8.74300632e-16, -9.07596236e-03, -1.85559638e-01],
                       [-1.09634524e-15, 5.89236266e-02, -3.46944487e-18],
                       [7.28583860e-17, 7.84809723e-02, 1.25453130e-01],
                       [7.97972799e-17, 7.84809723e-02, -1.25453130e-01],
                       [1.73472348e-16, -1.97733646e-01, -1.08422864e-18],
                       [-3.71924713e-15, -9.07596236e-03, 1.85559638e-01],
                       [-3.70536934e-15, -9.07596236e-03, -1.85559638e-01],
                       [1.14498795e-15, 5.89236266e-02, 6.93889390e-18],
                       [-2.08166817e-17, 7.84809723e-02, 1.25453130e-01],
                       [-2.08166817e-17, 7.84809723e-02, -1.25453130e-01]]
                      )

stress_ref = np.array([-2.23665496e+00, -1.44393236e+00, -1.97804121e+00,
                       5.55111512e-17, 2.77555756e-17, 4.71844785e-15])


modelname = 'SW_MX2_WenShirodkarPlechac_2017_MoS__MO_201919462778_001'
calc = KIM(modelname)

mos2 = ase.build.mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.18,
                     size=(2, 2, 1), vacuum=None)
mos2.set_pbc((1, 0, 0))

mos2.set_calculator(calc)

energy = mos2.get_potential_energy()
forces = mos2.get_forces()
stress = mos2.get_stress()

tol = 1e-6
assert np.isclose(energy, energy_ref, tol)
assert np.allclose(forces, forces_ref, tol)
assert np.allclose(stress, stress_ref, tol)
