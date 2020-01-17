import numpy as np
from ase.geometry.cell import cell_to_cellpar


def get_cell_angles_lengths(cell):
    '''
    Returns cell vectors lengths (a,b,c) as well as different
    angles (alpha, beta, gamma, phi, chi, psi) (in radians).
    '''
    cellpar = cell_to_cellpar(cell)
    cellpar[3:] *= np.pi / 180  # convert angles to radians
    parnames = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    values = {n: p for n, p in zip(parnames, cellpar)}

    volume = abs(np.linalg.det(cell))
    for i, param in enumerate(['phi', 'chi', 'psi']):
        ab = np.linalg.norm(
            np.cross(cell[(i + 1) % 3, :], cell[(i + 2) % 3, :]))
        c = np.linalg.norm(cell[i, :])
        s = np.abs(volume / (ab * c))
        if 1 + 1e-6 > s > 1:
            s = 1.
        values[param] = np.arcsin(s)

    return values


class CellBounds:
    '''
    Class for defining as well as checking limits on
    cell vector lengths and angles

    Parameters:

    bounds: dict
        Any of the following keywords can be used, in
        conjunction with a [low, high] list determining
        the lower and upper bounds:

        a, b, c:
           Minimal and maximal lengths (in Angstrom)
           for the 1st, 2nd and 3rd lattice vectors.

        alpha, beta, gamma:
           Minimal and maximal values (in degrees)
           for the angles between the lattice vectors.

        phi, chi, psi:
           Minimal and maximal values (in degrees)
           for the angles between each lattice vector
           and the plane defined by the other two vectors.

    Example:

    >>> from ase.ga.bulk_utilities import CellBounds
    >>> CellBounds(bounds={'phi': [20, 160],
    ...                    'chi': [60, 120],
    ...                    'psi': [20, 160],
    ...                    'a': [2, 20], 'b': [2, 20], 'c': [2, 20]})
    '''

    def __init__(self, bounds={}):
        self.bounds = {'alpha': [0, np.pi], 'beta': [0, np.pi],
                       'gamma': [0, np.pi], 'phi': [0, np.pi],
                       'chi': [0, np.pi], 'psi': [0, np.pi],
                       'a': [0, 1e6], 'b': [0, 1e6], 'c': [0, 1e6]}

        for param, bound in bounds.items():
            if param not in ['a', 'b', 'c']:
                # Convert angle from degree to radians
                bound = [b * np.pi / 180. for b in bound]
            self.bounds[param] = bound

    def is_within_bounds(self, cell):
        values = get_cell_angles_lengths(cell)
        verdict = True
        for param, bound in self.bounds.items():
            if not (bound[0] <= values[param] <= bound[1]):
                verdict = False
        return verdict


def get_rotation_matrix(u, t):
    '''
    Returns the transformation matrix for rotation over an angle t
    along an axis with direction u.
    '''
    ux, uy, uz = u
    cost, sint = np.cos(t), np.sin(t)
    rotmat = np.array([[(ux**2) * (1 - cost) + cost,
                        ux * uy * (1 - cost) - uz * sint,
                        ux * uz * (1 - cost) + uy * sint],
                       [ux * uy * (1 - cost) + uz * sint,
                        (uy**2) * (1 - cost) + cost,
                        uy * uz * (1 - cost) - ux * sint],
                       [ux * uz * (1 - cost) - uy * sint,
                        uy * uz * (1 - cost) + ux * sint,
                        (uz**2) * (1 - cost) + cost]])
    return rotmat
