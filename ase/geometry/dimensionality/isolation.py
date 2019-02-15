"""
Implements functions for extracting ('isolating') a low-dimensional material
component in its own unit cell.

Uses the rank-determination method described in:
Definition of a scoring parameter to identify low-dimensional materials
components
P.M. Larsen, M. Pandey, M. Strange, and K. W. Jacobsen
https://arxiv.org/abs/1808.02114
2018
"""


import itertools
import numpy as np

import ase
from ase.data import covalent_radii
from ase.neighborlist import NeighborList

from ase.geometry.dimensionality import interval_analysis
from ase.geometry.dimensionality import rank_determination
from ase.geometry.dimensionality.disjoint_set import DisjointSet


def orthogonal_basis(X, Y=None):

    is_1d = Y is None
    if Y is None:
        while 1:
            Y = np.random.uniform(-1, 1, 3)
            if abs(np.dot(X, Y)) < 0.5:
                break

    b = np.array([X, Y, np.cross(X, Y)])
    Q = np.linalg.qr(b.T)[0].T

    if np.dot(X, Q[0]) < 0:
        Q[0] = -Q[0]

    if np.dot(Y, Q[1]) < 0:
        Q[1] = -Q[1]

    if np.linalg.det(Q) < 0:
        Q[2] = -Q[2]

    if is_1d:
        Q = Q[[1, 2, 0]]

    return Q


def traverse_graph(atoms, k):

    num_atoms = len(atoms)
    rs = covalent_radii[atoms.get_atomic_numbers()]
    nl = NeighborList(k * rs, skin=0, self_interaction=False)
    nl.update(atoms)
    bonds = interval_analysis.get_bond_list(atoms, nl, rs)

    graph_bonds = []
    graph = DisjointSet(num_atoms)
    for (k, i, j, offset) in bonds:
        roffset = tuple(-np.array(offset))
        graph_bonds.append((i, j, offset))
        graph_bonds.append((j, i, roffset))
        if offset == (0, 0, 0):
            graph.merge(i, j)
    components = graph.get_components()

    adj = rank_determination.build_adjacency_list(components, graph_bonds)
    all_visited, ranks = rank_determination.traverse_component_graphs(adj)
    return components, all_visited


def build_supercomponent(atoms, components, k, v):

    # build superlayer/superchain by mapping components into visited cells
    positions = []
    numbers = []
    seen = set()
    for c, offset in v:

        # only want one copy of each sub-component
        if c in seen:
            continue
        else:
            seen.add(c)

        indices = np.where(components == c)[0]
        ps = atoms.positions[indices] + np.dot(offset, atoms.get_cell())
        positions += list(ps)
        numbers += list(atoms.numbers[indices])
    positions = np.array(positions)
    numbers = np.array(numbers)

    # select an 'anchor' atom, which will lie at the origin
    anchor_index = next((i for i in range(len(atoms)) if components[i] == k))
    anchor = atoms.positions[anchor_index]
    positions -= anchor
    return positions, numbers


def select_chain_rotation(scaled):

    best = (-1, None)
    for s in scaled:
        vhat = np.array([s[0], s[1], 0])
        vhat /= np.linalg.norm(vhat)
        obj = np.sum(np.dot(scaled, vhat)**2)
        best = max(best, (obj, vhat), key=lambda x: x[0])
    _, vhat = best
    cost, sint, _ = vhat
    rot = np.array([[cost, -sint, 0], [sint, cost, 0], [0, 0, 1]])
    return np.dot(scaled, rot)


def isolate_chain(atoms, components, k, v):

    positions, numbers = build_supercomponent(atoms, components, k, v)

    # identify the vector along the chain; this is the new cell vector
    basis_points = np.array([offset for c, offset in v if c == k])
    assert len(basis_points) >= 2
    assert (0, 0, 0) in [tuple(e) for e in basis_points]

    sizes = np.linalg.norm(basis_points, axis=1)
    index = np.argsort(sizes)[1]
    basis = basis_points[index]
    vector = np.dot(basis, atoms.get_cell())
    norm = np.linalg.norm(vector)
    vhat = vector / norm

    # project atoms into new basis
    scaled = np.dot(positions, orthogonal_basis(vhat).T / norm)

    # move atoms into new cell
    scaled[:,2] %= 1.0

    # subtract barycentre in x and y directions
    scaled[:, :2] -= np.mean(scaled, axis=0)[:2]

    # pick a good chain rotation (i.e. non-random) 
    scaled = select_chain_rotation(scaled)

    # construct a new atoms object containing the isolated chain
    cell = norm * np.eye(3)
    return ase.Atoms(numbers=numbers, scaled_positions=scaled, cell=cell,
                 pbc=[0, 0, 1])


def construct_inplane_basis(atoms, k, v):

    basis_points = np.array([offset for c, offset in v if c == k])
    assert len(basis_points) >= 3
    assert (0, 0, 0) in [tuple(e) for e in basis_points]

    sizes = np.linalg.norm(basis_points, axis=1)
    indices = np.argsort(sizes)
    basis_points = basis_points[indices]

    # identify primitive basis
    best = (float("inf"), None)
    for u, v in itertools.combinations(basis_points, 2):

        basis = np.array([[0, 0, 0], u, v])
        if np.linalg.matrix_rank(basis) < 2:
            continue

        a = np.dot(u, atoms.get_cell())
        b = np.dot(v, atoms.get_cell())
        norm = np.linalg.norm(np.cross(a, b))
        best = min(best, (norm, a, b), key=lambda x: x[0])
    _, a, b = best
    return a, b, orthogonal_basis(a, b)


def isolate_monolayer(atoms, components, k, v):

    positions, numbers = build_supercomponent(atoms, components, k, v)
    a, b, basis = construct_inplane_basis(atoms, k, v)

    # project atoms into new basis
    c = np.cross(a, b)
    c /= np.linalg.norm(c)
    cell = np.dot(np.array([a, b, c]), basis.T)
    scaled = np.linalg.solve(cell.T, np.dot(positions, basis.T).T).T

    # move atoms into new cell
    scaled[:,:2] %= 1.0

    # subtract barycentre in z direction
    scaled[:, 2] -= np.mean(scaled, axis=0)[2]

    # construct a new atoms object containing the isolated chain
    return ase.Atoms(numbers=numbers, scaled_positions=scaled,
                 cell=cell, pbc=[1, 1, 0])


def isolate_components(atoms, k):

    chains = {}
    monolayers = {}
    components, all_visited = traverse_graph(atoms, k)

    for k, v in all_visited.items():
        v = sorted(list(v))

        # identify the components which constitute the component
        key = tuple(np.unique([c for c, offset in v]))

        cells = np.array([offset for c, offset in v if c == k])
        rank = rank_determination.calc_rank(cells)

        if rank == 1:
            chains[key] = isolate_chain(atoms, components, k, v)
        elif rank == 2:
            monolayers[key] = isolate_monolayer(atoms, components, k, v)

    return list(chains.values()), list(monolayers.values())
