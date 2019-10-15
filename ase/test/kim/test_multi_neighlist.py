"""
Construct a 10 Angstrom x 10 Angstrom x 10 Angstrom non-periodic cell
filled with randomly positioned atoms and compute the energy, forces,
and virial stress for a model that makes use of multiple cutoffs
"""
import numpy as np
from ase import Atoms
from ase.calculators.kim import KIM

# Create random cluster of atoms
positions = np.random.RandomState(34).rand(15, 3) * 10
atoms = Atoms(
    "Ar" * 15, positions=positions, pbc=False, cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]]
)

calc = KIM("ex_model_Ar_P_Morse_MultiCutoff")
atoms.set_calculator(calc)

# Get energy and analytical forces/stress from KIM Model
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()

# Previously computed energy for this configuration for this model
energy_ref = 34.69963483186903

# Compute forces and virial stress numerically
forces_numer = calc.calculate_numerical_forces(atoms, d=0.0001)
stress_numer = calc.calculate_numerical_stress(atoms, d=0.0001, voigt=True)

tol = 1e-6
assert np.isclose(energy, energy_ref, tol)
assert np.allclose(forces, forces_numer, tol)
assert np.allclose(stress, stress_numer, tol)
