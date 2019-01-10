"""Implements the Topological Scaling Algorithm (TSA)

Method is described in:
Topology-Scaling Identification of Layered Solids and Stable Exfoliated
2D Materials
M. Ashton, J. Paul, S.B. Sinnott, and R.G. Hennig
Phys. Rev. Lett. 118, 106101
2017


A disjoint set is used here to allow insertion of bonds one at a time.
This permits k-interval analysis.
"""


import itertools
import numpy as np
from ase.geometry.dimensionality.disjoint_set import DisjointSet


class TSA:

    def __init__(self, num_atoms, n=2):

        """Initializes the TSA class.

        A disjoint set is maintained for the single cell and for the supercell.
        For some materials, such as interpenetrating networks,
        the dimensionality classification is dependent on the size of the
        initial cell.


        Parameters:

        num_atoms: int    The number of atoms in the unit cell.
        n: int            The number size of the (n, n, n) periodic supercell.
        """

        self.n = n
        self.num_atoms = num_atoms
        self.gsingle = DisjointSet(num_atoms)
        self.gsuper = DisjointSet(num_atoms * n**3)

        self.m = [1, n, n**2]
        self.cells = np.array(list(itertools.product(range(n), repeat=3)))
        self.offsets = num_atoms * np.dot(self.m, self.cells.T)

    def insert_bond(self, i, j, offset):

        """Inserts a bond into the component graph, both in the single cell and
        each of the n^3 subcells of the supercell.


        Parameters:

        i: int           The index of the first atom.
        n: int           The index of the second atom.
        offset: tuple    The cell offset of the second atom.
        """

        nbr_cells = (self.cells + offset) % self.n
        nbr_offsets = self.num_atoms * np.dot(self.m, nbr_cells.T)

        self.gsingle.merge(i, j)
        for (a, b) in zip(self.offsets, nbr_offsets):
            self.gsuper.merge(a + i, b + j)
            self.gsuper.merge(b + i, a + j)

    def check(self):

        """Determines the dimensionality of each component with the TSA method.

        Returns:

        hist : tuple         Dimensionality histogram.
        components: array    The component ID every atom
        """

        n = self.n
        offsets = self.offsets
        single_roots = self.gsingle.get_roots()
        super_components = self.gsuper.get_components()

        component_dim = {}
        hist = np.zeros(4).astype(np.int)
        for i in single_roots:

            num_clusters = len(np.unique(super_components[offsets + i]))
            dim = {n**3: 0, n**2: 1, n: 2, 1: 3}[num_clusters]
            hist[dim] += 1
            component_dim[i] = dim

        relabelled_components = self.gsingle.get_components(relabel=True)
        relabelled_dim = {}
        for k, v in component_dim.items():
            relabelled_dim[relabelled_components[k]] = v

        return tuple(hist), relabelled_components, relabelled_dim
