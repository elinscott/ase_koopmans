# See https://gitlab.com/ase/ase/issues/502
import numpy as np
from ase import Atoms
from ase.dft.kpoints import bandpath, labels_from_kpts

a = 3.16
atoms = Atoms(cell=[a, a, 12, 90, 90, 120], pbc=True)
G = [0, 0, 0]
K = [1 / 3., 1 / 3., 0]
K_ = [-1 / 3., -1 / 3., 0]
path = [K, G, K_]
nspecial_points = len(path)

for npoints in [10, 11]:
    bandpath = atoms.cell.bandpath(path, npoints=npoints)
    kpts, x, X = bandpath.get_linear_kpoint_axis()
    assert len(kpts) == npoints
    assert len(x) == nspecial_points
    assert len(X) == nspecial_points
