import re
from string import digits
import numpy as np
from ase import Atoms
from ase.units import Angstrom, Bohr, nm


_dist = {'angstrom': Angstrom, 'bohr': Bohr, 'au': Bohr, 'nm': nm}
_angle = {'radians': 1., 'degrees': np.pi / 180}


def _get_units(name, value, known_values):
    if isinstance(value, str):
        out = known_values.get(value.lower(), None)
        if out is None:
            raise ValueError("Unknown {} units: {}".format(name, value))
        return out
    else:
        return float(value)


# split on newlines or semicolons
_re_linesplit = re.compile(r'\n|;')
# split definitions on whitespace or on "=" (possibly also with whitespace)
_re_defs = re.compile(r'\s*=\s*|\s+')


def _get_var(key, defs, conv=1):
    # sometimes keys are given with a sign, e.g. "+L1" or "-A3".
    sign = -1 if key.startswith('-') else 1
    key = key.lstrip('+-')
    return sign * conv * float(defs.get(key, key))


def parse_zmatrix(zmat, distance_units='angstrom', angle_units='degrees',
                  defs=None):
    """Converts a Z-matrix into an Atoms object.

    Parameters:

    zmat: Iterable or str
        The Z-matrix to be parsed. Iteration over `zmat` should yield the rows
        of the Z-matrix. If `zmat` is a str, it will be automatically split
        into a list at newlines.
    distance_units: str or float, optional
        The units of distance in the provided Z-matrix.
        Defaults to Angstrom.
    angle_units: str or float, optional
        The units for angles in the provided Z-matrix.
        Defaults to degrees.
    defs: dict or str, optional
        If `zmat` contains symbols for bond distances, bending angles, and/or
        dihedral angles instead of numeric values, then the definition of
        those symbols should be passed to this function using this keyword
        argument.
        Note: The symbol definitions are typically printed adjacent to the
        Z-matrix itself, but this function will not automatically separate
        the symbol definitions from the Z-matrix.

    Returns:

    atoms: Atoms object
    """

    if defs is None:
        defs = dict()
    elif isinstance(defs, str):
        newdefs = dict()
        for row in _re_linesplit.split(defs.strip()):
            tokens = _re_defs.split(row.strip())
            if len(tokens) == 0:
                continue
            elif len(tokens) == 2:
                # don't convert to float here, that happens below
                newdefs[tokens[0]] = tokens[1]
            else:
                raise ValueError("Failed to parse the definition block! "
                                 "Problematic entry: {}".format(row))
        defs = newdefs

    dconv = _get_units('distance', distance_units, _dist)
    aconv = _get_units('angle', angle_units, _angle)

    # Start by assuming the z-matrix is 1-indexed, even though it probably
    # isn't. We figure that out below and modify offset apropriately.
    offset = 0

    # zmat should be a list containing the rows of the z-matrix.
    # for convenience, allow block strings and split at newlines.
    if isinstance(zmat, str):
        zmat = _re_linesplit.split(zmat.strip())

    natoms = len(zmat)
    positions = np.zeros((natoms, 3))

    # when the Z-matrix uses a distinct symbol for each atom (e..g 'C1', 'O2',
    # etc.), it is possible to refer to other atoms in the Z-matrix by their
    # symbol rather than their row number. In that case, we still need to know
    # the row number for the atoms, which is the purpose of this dict.
    atomname_to_index = dict()

    symbols = []
    for n, row in enumerate(zmat):
        tokens = row.strip().split()

        atomname_to_index[tokens[0]] = n

        # sometimes the element symbol has the atom index as well, e.g. "C1",
        # "O2", "O3". That's extra fun because "C1" could also refer to the
        # point group of the molecule. Let's just pretend that never happens.
        symbol = ''.join([char for char in tokens[0] if char not in digits])
        symbols.append(symbol.capitalize())

        if n == 0:
            continue

        a = atomname_to_index.get(tokens[1])
        if a is None:
            a = int(tokens[1]) + offset
        # Here and below, we verify the strict ordering of the atom indices
        # in the z-matrix. If we *don't* check, and the user provides an
        # invalid z-matrix, then either the user will get a completely bogus
        # structure, or the code further below will divide by zero, which
        # will raise error messages that will likely confuse users.
        if a >= n:
            if n == 1 and a == 1:
                # Turns out the z-matrix was 1-indexed. Let's fix offset.
                offset = -1
                a += offset
            else:
                raise ValueError("An invalid Z-matrix was provided!")
        dist = _get_var(tokens[2], defs, dconv)

        if n == 1:
            positions[n, 0] = dist
            continue

        b = atomname_to_index.get(tokens[3])
        if b is None:
            b = int(tokens[3]) + offset
        if b >= n or b == a:
            raise ValueError("An invalid Z-matrix was provided!")
        a_bend = _get_var(tokens[4], defs, aconv)

        if n == 2:
            vec = np.zeros(3)
            # this is a convenient trick: if n == 2, then either
            # (a, b) == (0, 1) OR (a, b) == (1, 0), so (b - a) == +/- 1
            # The 0->1 bond vector always points in the positive X direction,
            # so the a->n bond vector will start out pointing in the positive
            # X direction IFF a == 1, or in the negative X direction IFF
            # a == 0. (and then it will be rotated into the Y plane according
            # to the bending angle, a_bend).
            positions[n] = positions[a]
            positions[n, 0] += dist * np.cos(a_bend) * (b - a)
            positions[n, 1] += dist * np.sin(a_bend)
            continue

        c = atomname_to_index.get(tokens[5])
        if c is None:
            c = int(tokens[5]) + offset
        if c >= n or c == a or c == b:
            raise ValueError("An invalid Z-matrix was provided!")
        a_dihedral = _get_var(tokens[6], defs, aconv)

        # ax1 is the dihedral axis, which is defined by the bond vector
        # between the two inner atoms in the dihedral, a and b
        ax1 = positions[b] - positions[a]
        ax1 /= np.linalg.norm(ax1)

        # ax2 lies within the a-b-c plane, and it is perpendicular
        # to the dihedral axis
        ax2 = positions[b] - positions[c]
        ax2 -= ax1 * (ax2 @ ax1)
        ax2 /= np.linalg.norm(ax2)

        # ax3 is perpendicular to both ax2 and ax2.
        # The ax1, ax2, and ax3 form a complete basis in R^3
        ax3 = np.cross(ax2, ax1)
        ax3 /= np.linalg.norm(ax3)

        # ax4 is a vector that forms the appropriate dihedral angle, though
        # the bending angle is 90 degrees, rather than a_bend. It is formed
        # from a linear combination of ax2 and ax3.
        ax4 = ax2 * np.cos(a_dihedral) + ax3 * np.sin(a_dihedral)
        ax4 -= ax1 * (ax4 @ ax1)  # this step may not be necessary
        ax4 /= np.linalg.norm(ax4)

        # The final position vector is a linear combination of ax1 and ax4.
        vec = ax1 * np.cos(a_bend) - ax4 * np.sin(a_bend)
        vec *= dist / np.linalg.norm(vec)
        positions[n] = positions[a] + vec
    return Atoms(symbols, positions)
