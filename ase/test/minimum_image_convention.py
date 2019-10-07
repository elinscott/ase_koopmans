import numpy as np
from numpy.testing import assert_allclose
from ase.lattice.cubic import FaceCenteredCubic
from ase.neighborlist import mic as NeighborListMic


size = 2
atoms = FaceCenteredCubic(size=[size, size, size],
                          symbol='Cu',
                          latticeconstant=2,
                          pbc=(1, 1, 1))

d0 = atoms.get_distances(0, np.arange(len(atoms)), mic=True)

U = np.array([[1, 2, 2],
              [0, 1, 2],
              [0, 0, 1]])
assert np.linalg.det(U) == 1
atoms.set_cell(U.T @ atoms.cell, scale_atoms=False)
atoms.wrap()

d1 = atoms.get_distances(0, np.arange(len(atoms)), mic=True)
assert_allclose(d0, d1)

vnbrlist = NeighborListMic(atoms.get_positions(), atoms.cell, atoms.pbc)
d2 = np.linalg.norm(vnbrlist, axis=1)
assert_allclose(d0, d2)
