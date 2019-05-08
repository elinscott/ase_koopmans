"""
Implements the Rank Determination Algorithm (RDA)

Method is described in:
Definition of a scoring parameter to identify low-dimensional materials
components
P.M. Larsen, M. Pandey, M. Strange, and K. W. Jacobsen
Phys. Rev. Materials 3 034003, 2019
https://doi.org/10.1103/PhysRevMaterials.3.034003
"""

import numpy as np
from collections import defaultdict
from ase.geometry.dimensionality.disjoint_set import DisjointSet


def fcalc_rank(l):
    """Fast geometric rank calculation.

    Only for arrays whose elements are linearly independent.
    """
    return len(l) - 1


def calc_rank(l):
    """Full rank calculation.

    The rank of an empty set of vectors is defined as -1.  The geometric rank
    of single point is zero.  The geometric rank of two linearly independent
    points (a line) is 1 etc.
    """
    return np.linalg.matrix_rank(np.array(l) - l[0])


def bfs(adjacency, start):
    """Traverse the component graph using BFS.

    The graph is traversed until the matrix rank of the subspace spanned by
    the visited components no longer increases.
    """

    visited = set()
    cvisited = defaultdict(list)
    queue = [(start, (0, 0, 0))]
    while queue:
        vertex = queue.pop(0)
        if vertex in visited:
            continue

        visited.add(vertex)
        c, pos = vertex
        if calc_rank(cvisited[c] + [pos]) > fcalc_rank(cvisited[c]):
            cvisited[c] += [pos]

        for nc, offset in adjacency[c]:

            nbrpos = tuple(np.array(pos) + offset)
            nbrnode = (nc, nbrpos)
            if nbrnode in visited:
                continue

            rank0 = fcalc_rank(cvisited[nc])
            rank1 = calc_rank(cvisited[nc] + [nbrpos])
            if rank1 == rank0:
                continue

            queue.append(nbrnode)

    return visited, fcalc_rank(cvisited[start])


def traverse_component_graphs(adjacency):

    vertices = adjacency.keys()
    all_visited = {}
    ranks = {}
    for v in vertices:
        visited, rank = bfs(adjacency, v)
        all_visited[v] = visited
        ranks[v] = rank

    return all_visited, ranks


def build_adjacency_list(parents, bonds):

    graph = np.unique(parents)
    adjacency = {e: set() for e in graph}
    for (i, j, offset) in bonds:
        component_a = parents[i]
        component_b = parents[j]
        adjacency[component_a].add((component_b, offset))
    return adjacency


def get_dimensionality_histogram(ranks, roots):

        component_ranks = [ranks[e] for e in roots]

        h = np.zeros(4).astype(np.int)
        bc = np.bincount(component_ranks)
        h[:len(bc)] = bc
        return tuple(h)


class RDA:

    def __init__(self, num_atoms):

        """
        Initializes the RDA class.

        A disjoint set is used to maintain the component graph.

        Parameters:

        num_atoms: int    The number of atoms in the unit cell.
        """

        self.bonds = []
        self.graph = DisjointSet(num_atoms)
        self.adjacency = None
        self.hcached = None
        self.components_cached = None
        self.cdim_cached = None

    def insert_bond(self, i, j, offset):

        """
        Adds a bond to the list of graph edges.

        Graph components are merged if the bond does not cross a cell boundary.
        Bonds which cross cell boundaries can inappropriately connect
        components which are not connected in the infinite crystal.  This is
        tested during graph traversal.
        

        Parameters:

        i: int           The index of the first atom.
        n: int           The index of the second atom.
        offset: tuple    The cell offset of the second atom.
        """

        roffset = tuple(-np.array(offset))
        self.bonds += [(i, j, offset)]
        self.bonds += [(j, i, roffset)]

        if offset == (0, 0, 0):    # only want bonds in aperiodic unit cell
            self.graph.merge(i, j)

    def check(self):

        """
        Determines the dimensionality of each component using the RDA method.

        The component graph is traversed (using BFS) until the matrix rank
        of the subspace spanned by the visited components no longer increases.

        Returns:
        hist : tuple         Dimensionality histogram.
        components: array    The component ID every atom
        """

        adjacency = build_adjacency_list(self.graph.get_components(),
                                         self.bonds)
        if adjacency == self.adjacency:
            return self.hcached, self.components_cached, self.cdim_cached

        self.adjacency = adjacency
        all_visited, ranks = traverse_component_graphs(adjacency)

        # Find components with mutual visits
        common = defaultdict(list)
        for c, visited in all_visited.items():
            for e in visited:
                common[e] += [c]

        # Merge components with mutual visits
        for k, v in common.items():
            for i in range(len(v) - 1):
                a = v[i]
                b = v[i + 1]
                assert ranks[a] == ranks[b]
                self.graph.merge(a, b)

        h = get_dimensionality_histogram(ranks, self.graph.get_roots())
        self.hcached = h

        component_dim = {e: ranks[e] for e in self.graph.get_roots()}
        relabelled_components = self.graph.get_components(relabel=True)
        relabelled_dim = {}
        for k, v in component_dim.items():
            relabelled_dim[relabelled_components[k]] = v
        self.cdim_cached = relabelled_dim
        self.components_cached = relabelled_components

        return h, relabelled_components, relabelled_dim
