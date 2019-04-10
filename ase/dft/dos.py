import functools
from math import pi, sqrt

import numpy as np

from ase.dft.kpoints import get_monkhorst_pack_size_and_offset
from ase.parallel import world
from ase.utils.cext import cextension


class DOS:
    def __init__(self, calc, width=0.1, window=None, npts=401):
        """Electronic Density Of States object.

        calc: calculator object
            Any ASE compliant calculator object.
        width: float
            Width of guassian smearing.  Use width=0.0 for linear tetrahedron
            interpolation.
        window: tuple of two float
            Use ``window=(emin, emax)``.  If not specified, a window
            big enough to hold all the eigenvalues will be used.
        npts: int
            Number of points.

        """

        self.npts = npts
        self.width = width
        self.w_k = calc.get_k_point_weights()
        self.nspins = calc.get_number_of_spins()
        self.e_skn = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                                for k in range(len(self.w_k))]
                               for s in range(self.nspins)])
        self.e_skn -= calc.get_fermi_level()

        if window is None:
            emin = None
            emax = None
        else:
            emin, emax = window

        if emin is None:
            emin = self.e_skn.min() - 5 * self.width
        if emax is None:
            emax = self.e_skn.max() + 5 * self.width

        self.energies = np.linspace(emin, emax, npts)

        if width == 0.0:
            bzkpts = calc.get_bz_k_points()
            size, offset = get_monkhorst_pack_size_and_offset(bzkpts)
            bz2ibz = calc.get_bz_to_ibz_map()
            shape = (self.nspins,) + tuple(size) + (-1,)
            self.e_skn = self.e_skn[:, bz2ibz].reshape(shape)
            self.cell = calc.atoms.cell

    def get_energies(self):
        """Return the array of energies used to sample the DOS.

        The energies are reported relative to the Fermi level.
        """
        return self.energies

    def delta(self, energy):
        """Return a delta-function centered at 'energy'."""
        x = -((self.energies - energy) / self.width)**2
        return np.exp(x) / (sqrt(pi) * self.width)

    def get_dos(self, spin=None):
        """Get array of DOS values.

        The *spin* argument can be 0 or 1 (spin up or down) - if not
        specified, the total DOS is returned.
        """

        if spin is None:
            if self.nspins == 2:
                # Spin-polarized calculation, but no spin specified -
                # return the total DOS:
                return self.get_dos(spin=0) + self.get_dos(spin=1)
            else:
                spin = 0

        if self.width == 0.0:
            dos = linear_tetrahedron_integration(self.cell, self.e_skn[spin],
                                                 self.energies)
            return dos[0]

        dos = np.zeros(self.npts)
        for w, e_n in zip(self.w_k, self.e_skn[spin]):
            for e in e_n:
                dos += w * self.delta(e)
        return dos


def linear_tetrahedron_integration(cell, eigs, energies, weights=None):
    """DOS from linear tetrahedron interpolation.

    cell: 3x3 ndarray-like
        Unit cell.
    eigs: (n1, n2, n3, nbands)-shaped ndarray
        Eigenvalues on a Monkhorst-Pack grid (not reduced).
    energies: 1-d array-like
        Energies where the DOS is calculated (must be a uniform grid).
    weights: (n1, n2, n3, nbands, nweights)-shaped ndarray
        Weights.  Defaults to a (n1, n2, n3, nbands, 1)-shaped ndarray
        filled with ones.
    """

    from scipy.spatial import Delaunay

    size = eigs.shape[:3]
    B = (np.linalg.inv(cell) / size).T
    indices = np.array([[i, j, k]
                        for i in [0, 1] for j in [0, 1] for k in [0, 1]])
    dt = Delaunay(np.dot(indices, B))

    if weights is None:
        weights = np.ones_like(eigs)[..., np.newaxis]

    nweights = len(weights)
    dos = np.zeros((nweights, len(energies)))
    ltidos(indices[dt.simplices], eigs, weights, energies, dos, world)
    return dos


@cextension
def ltidos(simplices, eigs, weights, energies, dos, world):
    I, J, K = eigs.shape[:3]
    n = -1
    for i in range(I):
        for j in range(J):
            for k in range(K):
                for indices in simplices:
                    print(indices)
                    n += 1
                    if n % world.size != world.rank:
                        continue

                    E = []
                    W = []
                    for a, b, c in indices:
                        E.append(eigs[(i + a) % I,
                                      (j + b) % J,
                                      (k + c) % K])
                        W.append(weights[(i + a) % I,
                                         (j + b) % J,
                                         (k + c) % K])
                    for e1, e2, e3, e4, w1, w2, w3, w4 in zip(*E, *W):
                        ltidos1(sorted([(e1, w1), (e2, w2),
                                        (e3, w3), (e4, w4)]),
                                energies, dos)
    world.sum(dos)


def ltidos1(ew, energies, dos):
    print(ew)

