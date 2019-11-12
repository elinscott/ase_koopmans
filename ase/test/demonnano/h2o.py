from ase.calculators.demonnano import DemonNano
from ase import Atoms
from ase.optimize import BFGS
from ase.units import Bohr, Hartree
import numpy as np

d = 0.9775
t = np.pi / 180 * 110.51
atoms = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)])

input_arguments = {'DFTB': 'SCC',
                   'CHARGE': '0.0',
                   'PARAM': 'PTYPE=MAT'}
    
calc = DemonNano(input_arguments=input_arguments)

atoms.set_calculator(calc)

# energy
energy = atoms.get_potential_energy()
ref = -4.08209

print('energy')
print(energy)

error = np.sqrt(np.sum((energy - ref)**2))
print('diff from reference:')
print(error)

tol = 1.0e-6
assert(error < tol)

# analytical forces
forces_an = atoms.get_forces()
ref = np.array([[ 0.1123E-01, -0.1670E-01, 0.0000E+00],
                [-0.1958E-01,  0.4671E-02, 0.0000E+00],
                [ 0.8349E-02,  0.1203E-01, 0.0000E+00]])
ref*=-Hartree/Bohr

error = np.sqrt(np.sum((forces_an - ref)**2))
print('forces_an')
print(forces_an)
print('diff from reference:')
print(error)

tol = 1.0e-3
assert(error < tol)

# optimize geometry
dyn = BFGS(atoms)
dyn.run(fmax=0.01)

positions = atoms.get_positions()

ref = np.array([[ 0.943765,   0.046188,   0.000000],
                [-0.287409,   0.900126,   0.000000],
                [-0.021346,  -0.030774,   0.000000]])

error = np.sqrt(np.sum((positions - ref)**2))
print('positions')
print(positions)
print('diff from reference:')
print(error)

tol = 1.0e-3
assert(error < tol)

print('tests passed')
