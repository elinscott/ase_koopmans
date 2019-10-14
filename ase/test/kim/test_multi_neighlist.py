"""
Construct a 10 Angstrom x 10 Angstrom x 10 Angstrom non-periodic cell
filled with randomly positioned atoms and compute the energy, forces,
and virial stress for a model that makes use of multiple cutoffs
"""
import numpy as np
from ase import Atoms
from ase.calculators.kim import KIM

# Set up a cluster of atoms
positions = np.random.RandomState(34).rand(15, 3) * 10

energy_ref = 34.69963483186903

forces_ref = np.array(
    [
        [-2.39112996e-02, 7.02682418e-02, -2.53387956e-02],
        [-4.20585468e-01, -2.58702688e-01, 1.95473299e-01],
        [-2.18260158e01, -6.29033132e01, 5.11156248e01],
        [9.16554524e00, -8.72108629e00, -4.72110931e00],
        [5.43675529e-02, -1.13377152e-01, 4.09716280e-02],
        [-1.73194441e00, 3.02473249e00, -1.34370347e00],
        [-9.25117418e00, 8.70590242e00, 4.59994850e00],
        [8.56127093e-02, 2.21861418e-01, 2.46975134e-01],
        [2.21959249e00, -2.96531047e00, 1.26753960e00],
        [-4.71392294e-01, 1.21182485e-01, 8.49547580e-02],
        [-2.38197176e-03, 2.30955640e-02, 8.65437764e-03],
        [9.96945240e-02, -5.35720199e-02, -3.08467397e-01],
        [-1.32884747e-01, 3.58642822e-01, 2.48936487e-01],
        [2.22392422e01, 6.24872421e01, -5.14150943e01],
        [-3.76454535e-03, 2.43426593e-03, 4.63467673e-03],
    ]
)

stress_ref = np.array(
    [-0.02526133, -0.08248774, -0.04771512, 0.04723221, 0.02290316, -0.00677126]
)

modelname = "ex_model_Ar_P_Morse_MultiCutoff"

calc = KIM(modelname)

atoms = Atoms(
    "Ar" * 15, positions=positions, pbc=False, cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]]
)
atoms.set_calculator(calc)

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()

tol = 1e-6
assert np.isclose(energy, energy_ref, tol)
assert np.allclose(forces, forces_ref, tol)
assert np.allclose(stress, stress_ref, tol)
