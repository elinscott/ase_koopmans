"""Implements the dimensionality scoring parameter.

Method is described in:
Definition of a scoring parameter to identify low-dimensional materials
components
P.M. Larsen, M. Pandey, M. Strange, and K. W. Jacobsen
https://arxiv.org/abs/1808.02114
2018
"""

import numpy as np
from ase.neighborlist import NeighborList
from ase.data import covalent_radii as radii
from collections import defaultdict
from ase.geometry.dimensionality import rank_determination
from ase.geometry.dimensionality import topological_scaling


def f(x):
    if x == float("inf"):
        return 1
    k = 1 / 0.15**2
    return k * max(0, x - 1)**2 / (1. + k * max(0, x - 1)**2)


def calculate_score(a, b):
    return f(b) - f(a)


def reduced_histogram(h):

    h = [int(e > 0) for e in h]
    return tuple(h)


def hstring(h):

    h = reduced_histogram(h)
    return ''.join([str(i) for i, e in enumerate(h) if e > 0]) + 'D'


def merge_intervals(intervals):

    """Merges intervals of the same type (same reduced histogram).

    For example, two histograms with component histograms [10, 4, 0, 0] and
    [6, 2, 0, 0] are both 01D structures so they will be merged.

    Intervals are merged by summing the scores, and taking the minimum and
    maximum k-values.  Component IDs in the merged interval are taken from the
    interval with the highest score.

    On rare occasions, intervals to be merged are not adjacent.  In this case,
    the score of the merged interval is not equal to the score which would be
    calculated from its k-interval.  This is necessary to maintain the property
    that the scores sum to 1.
    """

    merged = defaultdict(list)
    for (a, b, h, components, cdim) in intervals:
        hr = hstring(h)
        merged[hr] += [(b - a, a, b, components, h, cdim)]

    merged_intervals = []
    for hr, intervals in merged.items():
        amin = min([a for _, a, b, _, _, _ in intervals])
        bmax = max([b for _, a, b, _, _, _ in intervals])
        score = sum([calculate_score(a, b) for _, a, b, _, _, _ in intervals])
        _, _, _, components, h, cdim = max(intervals)
        merged_intervals += [(score, amin, bmax, hr, h, components, cdim)]

    return sorted(merged_intervals, reverse=True)


def get_bond_list(atoms, nl, rs):

    """Gets a list of bonds sorted by k-value, from low to high.

    Parameters:

    atoms: ASE atoms object
        The periodic solid to analyze
    nl: ASE neighborlist
    rs: Scaled covalent radii (k-value times covalent radii)

    Returns:

    intervals : list
        List of tuples for each bond.  Each tuple contains
        (k, i, j, offset)

        k:       float   k-value
        i:       float   index of first atom
        j:       float   index of second atom
        offset:  tuple   cell offset of second atom
    """

    num_atoms = len(atoms)
    bonds = []
    for i in range(num_atoms):
        p = atoms.positions[i]
        indices, offsets = nl.get_neighbors(i)

        for j, offset in zip(indices, offsets):
            q = atoms.positions[j] + np.dot(offset, atoms.get_cell())
            d = np.linalg.norm(p - q)
            k = d / (rs[i] + rs[j])
            bonds += [(k, i, j, tuple(offset))]
    return sorted(bonds)


def analyze_kintervals(atoms, method='RDA'):

    """Performs a k-interval analysis of a periodic solid.

    In each k-interval the components (connected clusters) are identified.
    The intervals are sorted according to the scoring parameter, from high
    to low.


    Parameters:

    atoms: ASE atoms object
        The periodic solid to analyze
    method: string
        Analysis method to use, either 'RDA' (default option) or 'TSA'.
        These correspond to the Rank Determination Algorithm of Mounet et al.
        and the Topological Scaling Algorithm (TSA) of Ashton et al.

    Returns:

    intervals : list
        List of tuples for each interval identified.  Each tuple contains
        (score, amin, bmax, hr, components)

        score:      float   Dimensionality score in the range [0, 1].
                            A higher score is better.
        amin:       float   The start of the k-interval
        bmax:       float   The end of the k-interval
        hr:         string  The reduced histogram of the interval
                            e.g. 0D, 1D, 2D, 3D, 03D, 012D etc.
        h:          tuple   The histogram of the number of components
                            of each dimensionality type, e.g. (0, 3, 0, 0)
                            means three 1D components.
        components: array   The component ID of each atom.
        cdim:       dict    The component dimensionalities
    """

    if method == 'RDA':
        method = rank_determination.RDA
    elif method == 'TSA':
        method = topological_scaling.TSA

    # Interval analysis requires a periodic solid
    assert tuple(atoms.pbc) == (1, 1, 1)
    num_atoms = len(atoms)
    rs = radii[atoms.get_atomic_numbers()]

    """
    The interval analysis is performed by iteratively expanding the neighbor
    lists, until the component analysis finds a single 3D component.  To avoid
    repeat analyses after expanding the neighbor lists, we keep track of the
    previously inserted bonds.
    """

    intervals = []
    seen = set()
    kprev = 0
    calc = method(num_atoms)
    hprev, components_prev, cdim_prev = calc.check()

    kmax = 0
    while 1:

        # Expand the scope of the neighbor lists.
        kmax += 2
        nl = NeighborList(kmax * rs, skin=0, self_interaction=False)
        nl.update(atoms)

        # Get a list of bonds, sorted by k-value.
        bonds = get_bond_list(atoms, nl, rs)

        # Find only the bonds which we have not previously tested.
        new_bonds = []
        for b in bonds:
            if b not in seen:
                new_bonds += [b]
                seen.add(b)

        # Insert each new bond into the component graph.
        for (k, i, j, offset) in new_bonds:

            calc.insert_bond(i, j, offset)
            h, components, cdim = calc.check()
            if h == hprev:    # Test if any components were merged
                continue

            # If any components were merged, create a new interval
            intervals += [(kprev, k, hprev, tuple(components_prev), cdim_prev)]
            kprev = k
            hprev = h
            components_prev = components
            cdim_prev = cdim

            # Stop once all components are merged (a single 3D component)
            if h == (0, 0, 0, 1):
                intervals += [(k, float("inf"), h, tuple(components), cdim)]
                return merge_intervals(intervals)
