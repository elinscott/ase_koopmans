import numpy as np
from ase.build import bulk, make_supercell

a = 4.0
atoms = bulk('Au', a=a)
assert atoms.cell.get_bravais_lattice().name == 'FCC'

P = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
assert abs(np.linalg.det(P) - 4) < 1e-14

cubatoms = make_supercell(atoms, P)
assert cubatoms.cell.orthorhombic
assert len(cubatoms) == 4
assert np.allclose(cubatoms.cell.lengths(), a)
assert cubatoms.cell.get_bravais_lattice().name == 'CUB'
