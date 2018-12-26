from __future__ import print_function, division
# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

# XXX bravais objects need to hold tolerance eps, *or* temember variant
# from the beginning.
#
# Should they hold a 'cycle' argument or other data to reconstruct a particular
# cell?  (E.g. rotation, niggli transform)
#
# Implement total ordering of Bravais classes 1-14

from abc import abstractmethod, ABC
import os
import re
import numpy as np
from numpy import pi, sin, cos, arccos, sqrt, dot
from numpy.linalg import norm

from ase.utils.arraywrapper import arraylike


@arraylike
class Cell:
    """Parallel epipedal unit cell of up to three dimensions.

    This wraps a 3x3 array whose [i, j]-th element is the jth
    Cartesian coordinate of the ith unit vector.

    Cells of less than three dimensions are represented by placeholder
    unit vectors that are zero."""

    # This overridable variable tells an Atoms object whether atoms.cell
    # and atoms.get_cell() should be a Cell object or an array.
    _atoms_use_cellobj = 1#bool(os.environ.get('ASE_DEBUG_CELLOBJ'))

    def __init__(self, array=None, pbc=None):
        if array is None:
            array = np.zeros((3, 3))

        if pbc is None:
            pbc = np.ones(3, bool)

        # We could have lazy attributes for structure (bcc, fcc, ...)
        # and other things.  However this requires making the cell
        # array readonly, else people will modify it and things will
        # be out of synch.
        assert array.shape == (3, 3)
        assert array.dtype == float
        assert pbc.shape == (3,)
        assert pbc.dtype == bool
        self.array = array
        self.pbc = pbc

    def cellpar(self, radians=False):
        return cell_to_cellpar(self.array, radians)

    @property
    def shape(self):
        return self.array.shape

    @classmethod
    def new(cls, cell):
        cell = np.array(cell, float)

        if cell.shape == (3,):
            cell = np.diag(cell)
        elif cell.shape == (6,):
            cell = cellpar_to_cell(cell)
        elif cell.shape != (3, 3):
            raise ValueError('Cell must be length 3 sequence, length 6 '
                             'sequence or 3x3 matrix!')

        return cls(cell)

    @classmethod
    def fromcellpar(cls, cellpar, ab_normal=(0, 0, 1), a_direction=None):
        cell = cellpar_to_cell(cellpar, ab_normal, a_direction)
        return Cell(cell)

    #def crystal_structure(self, eps=2e-4, niggli_reduce=True):
    #    return crystal_structure_from_cell(self.array, eps, niggli_reduce)

    def bravais(self, eps=2e-4):
        return get_bravais_lattice(self, eps=eps)

    def complete(self):
        """Convert missing cell vectors into orthogonal unit vectors."""
        return Cell(complete_cell(self.array))

    def copy(self):
        return Cell(self.array.copy())

    @property
    def dtype(self):
        return self.array.dtype

    @property
    def size(self):
        return self.array.size

    @property
    def T(self):
        return self.array.T

    @property
    def flat(self):
        return self.array.flat

    @property
    def celldim(self):
        # XXX Would name it ndim, but this clashes with ndarray.ndim
        return self.array.any(1).sum()

    @property
    def orthorhombic(self):
        return orthorhombic(self.array)

    @property
    def is_orthorhombic(self):
        return is_orthorhombic(self.array)

    @property
    def ndim(self):
        return self.array.ndim

    def box(self):
        """Return cell lengths if orthorhombic, else raise ValueError."""
        # XXX More intelligent name for thos method?
        return orthorhombic(self.array)

    def __array__(self, dtype=float):
        if dtype != float:
            raise ValueError('Cannot convert cell to array of type {}'
                             .format(dtype))
        return self.array

    def __bool__(self):
        return bool(self.array.any())

    def __ne__(self, other):
        return self.array != other

    def __eq__(self, other):
        return self.array == other

    __nonzero__ = __bool__

    @property
    def volume(self):
        # Fail or 0 for <3D cells?
        # Definitely 0 since this is currently a property.
        # I think normally it is more convenient just to get zero
        return np.abs(np.linalg.det(self.array))

    def tolist(self):
        return self.array.tolist()

    def scaled_positions(self, positions):
        return np.linalg.solve(self.complete().array.T, positions.T).T

    def cartesian_positions(self, scaled_positions):
        return np.dot(scaled_positions, self.complete().array)

    def reciprocal(self):
        return np.linalg.pinv(self.array).transpose()

    #def lattice_type(self, eps=2e-4):
    #    from ase.geometry.crystal_info import analyse_cell
    #    return analyse_cell(self, eps=eps)

    def __repr__(self):
        if self.is_orthorhombic:
            numbers = self.box().tolist()
        else:
            numbers = self.tolist()

        pbc = self.pbc
        if all(pbc):
            pbc = True
        elif not any(pbc):
            pbc = False
        return 'Cell({}, pbc={})'.format(numbers, pbc)

    def niggli_reduce(self):
        from ase.build.tools import niggli_reduce_cell
        cell, _ = niggli_reduce_cell(self.array)
        return Cell(cell)

    #def bandpath(self, path, npoints=50):
    #    from ase.dft.kpoints import bandpath, BandPath
    #    objs = bandpath(path, self.array, npoints=npoints)
    #    return BandPath(*objs, names=path)

    #def special_points(self, eps=2e-4):
    #    from ase.dft.kpoints import get_special_points
    #    return get_special_points(self.array, eps=eps)

    #def special_paths(self, eps=2e-4):
    #    from ase.dft.kpoints import special_paths
    #    structure = self.crystal_structure(eps=eps)
    #    pathstring = special_paths[structure]
    #    paths = pathstring.split(',')
    #    return paths

    #def bravais_type(self, eps=2e-4):
    #    """Get bravais lattice and lattice parameters as (lattice, par).

    #    lattice, par = uc.bravais()
    #    print(uc.cellpar())
    #    print(lattice(**par).cellpar())"""
    #    return get_bravais_lattice(self, eps=eps)


def unit_vector(x):
    """Return a unit vector in the same direction as x."""
    y = np.array(x, dtype='float')
    return y / norm(y)


def angle(x, y):
    """Return the angle between vectors a and b in degrees."""
    return arccos(dot(x, y) / (norm(x) * norm(y))) * 180. / pi


def cell_to_cellpar(cell, radians=False):
    """Returns the cell parameters [a, b, c, alpha, beta, gamma].

    Angles are in degrees unless radian=True is used.
    """
    lengths = [np.linalg.norm(v) for v in cell]
    angles = []
    for i in range(3):
        j = i - 1
        k = i - 2
        ll = lengths[j] * lengths[k]
        if ll > 1e-16:
            x = np.dot(cell[j], cell[k]) / ll
            angle = 180.0 / pi * arccos(x)
        else:
            angle = 90.0
        angles.append(angle)
    if radians:
        angles = [angle * pi / 180 for angle in angles]
    return np.array(lengths + angles)


def cellpar_to_cell(cellpar, ab_normal=(0, 0, 1), a_direction=None):
    """Return a 3x3 cell matrix from cellpar=[a,b,c,alpha,beta,gamma].

    Angles must be in degrees.

    The returned cell is orientated such that a and b
    are normal to `ab_normal` and a is parallel to the projection of
    `a_direction` in the a-b plane.

    Default `a_direction` is (1,0,0), unless this is parallel to
    `ab_normal`, in which case default `a_direction` is (0,0,1).

    The returned cell has the vectors va, vb and vc along the rows. The
    cell will be oriented such that va and vb are normal to `ab_normal`
    and va will be along the projection of `a_direction` onto the a-b
    plane.

    Example:

    >>> cell = cellpar_to_cell([1, 2, 4, 10, 20, 30], (0, 1, 1), (1, 2, 3))
    >>> np.round(cell, 3)
    array([[ 0.816, -0.408,  0.408],
           [ 1.992, -0.13 ,  0.13 ],
           [ 3.859, -0.745,  0.745]])

    """
    if a_direction is None:
        if np.linalg.norm(np.cross(ab_normal, (1, 0, 0))) < 1e-5:
            a_direction = (0, 0, 1)
        else:
            a_direction = (1, 0, 0)

    # Define rotated X,Y,Z-system, with Z along ab_normal and X along
    # the projection of a_direction onto the normal plane of Z.
    ad = np.array(a_direction)
    Z = unit_vector(ab_normal)
    X = unit_vector(ad - dot(ad, Z) * Z)
    Y = np.cross(Z, X)

    # Express va, vb and vc in the X,Y,Z-system
    alpha, beta, gamma = 90., 90., 90.
    if isinstance(cellpar, (int, float)):
        a = b = c = cellpar
    elif len(cellpar) == 1:
        a = b = c = cellpar[0]
    elif len(cellpar) == 3:
        a, b, c = cellpar
    else:
        a, b, c, alpha, beta, gamma = cellpar

    # Handle orthorhombic cells separately to avoid rounding errors
    eps = 2 * np.spacing(90.0, dtype=np.float64)  # around 1.4e-14
    # alpha
    if abs(abs(alpha) - 90) < eps:
        cos_alpha = 0.0
    else:
        cos_alpha = cos(alpha * pi / 180.0)
    # beta
    if abs(abs(beta) - 90) < eps:
        cos_beta = 0.0
    else:
        cos_beta = cos(beta * pi / 180.0)
    # gamma
    if abs(gamma - 90) < eps:
        cos_gamma = 0.0
        sin_gamma = 1.0
    elif abs(gamma + 90) < eps:
        cos_gamma = 0.0
        sin_gamma = -1.0
    else:
        cos_gamma = cos(gamma * pi / 180.0)
        sin_gamma = sin(gamma * pi / 180.0)

    # Build the cell vectors
    va = a * np.array([1, 0, 0])
    vb = b * np.array([cos_gamma, sin_gamma, 0])
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz = sqrt(1. - cx * cx - cy * cy)
    vc = c * np.array([cx, cy, cz])

    # Convert to the Cartesian x,y,z-system
    abc = np.vstack((va, vb, vc))
    T = np.vstack((X, Y, Z))
    cell = dot(abc, T)

    return cell


def metric_from_cell(cell):
    """Calculates the metric matrix from cell, which is given in the
    Cartesian system."""
    cell = np.asarray(cell, dtype=float)
    return np.dot(cell, cell.T)


class BravaisLattice(ABC):
    # These parameters can be set by the @bravais decorator for a subclass.
    # (We could also use metaclasses to do this, but that's more abstract)
    type = None  # e.g. 'CUB', 'BCT', 'ORCF', ...
    name = None  # e.g. 'cubic', 'body-centered tetragonal', ...
    parameters = None  # e.g. ('a', 'c')
    variants = None  # e.g. {'BCT1': <variant object>,
    #                        'BCT2': <variant object>}

    _class_desc_template = """\
Lattice name: {type}
  Long name: {name}
  Parameters: {parameters}
"""

    def __init__(self, **kwargs):
        p = {}
        for k, v in kwargs.items():
            p[k] = float(v)
        self._parameters = p

    def __getattr__(self, name):
        return self._parameters[name]

    def tocell(self, cycle=0):
        cell = self._cell(**self._parameters)
        if cycle:
            index = (np.arange(3) + cycle) % 3
            cell = cell[index]
        return Cell(cell)

    def cellpar(self, cycle=0):
        # (Just a brute-force implementation)
        cell = self.tocell(cycle=cycle)
        return cell.cellpar()

    def get_special_points(self):
        kw = self._parameters
        variant_name = self._variant_name(**kw)
        variant = self.variants[variant_name]
        points = self._special_points(variant=variant, **kw)
        # replace the string by list of points using some regex
        #assert len(points) == len(variant.special_point_names)
        return np.array(points)

    def get_variant(self):
        name = self._variant_name(**self._parameters)
        return self.variants[name]

    def get_kpoint_labels(self):
        variant = self.get_variant()
        labels = re.findall(r'[A-Z]\d?', variant.special_point_names)
        return labels

    @abstractmethod
    def _cell(self, **kwargs):
        pass

    @abstractmethod
    def _special_points(self, **kwargs):
        pass

    @abstractmethod
    def _variant_name(self, **kwargs):
        pass

    def __repr__(self):
        par = ', '.join('{}={}'.format(k, v)
                        for k, v in self._parameters.items())
        return '{}({})'.format(self.type, par)

    def __str__(self):
        special_points = self.get_special_points()

        labels = self.get_kpoint_labels()
        assert len(labels) == len(special_points)

        coordstring = '\n'.join(['    {:2s} {:7.4f} {:7.4f} {:7.4f}'
                                 .format(label, *point)
                                 for label, point
                                 in zip(labels, special_points)])

        string = """\
{repr}
  {variant}
  Special point coordinates:
{special_points}
""".format(repr=repr(self),
           variant=self.get_variant(),
           special_points=coordstring)
        return string

    @classmethod
    def type_description(cls):
        desc = cls._class_desc_template.format(**vars(cls))

        chunks = [desc]
        for name in cls.variant_names:
            var = cls.variants[name]
            txt = str(var)
            lines = ['  ' + L for L in txt.splitlines()]
            lines.append('')
            chunks.extend(lines)

        return '\n'.join(chunks)


ibz_points = {'cub': {'G': [0, 0, 0],
                      'X': [0, 0 / 2, 1 / 2],
                      'R': [1 / 2, 1 / 2, 1 / 2],
                      'M': [0 / 2, 1 / 2, 1 / 2]},
              'fcc': {'G': [0, 0, 0],
                      'X': [1 / 2, 0, 1 / 2],
                      'W': [1 / 2, 1 / 4, 3 / 4],
                      'K': [3 / 8, 3 / 8, 3 / 4],
                      'U': [5 / 8, 1 / 4, 5 / 8],
                      'L': [1 / 2, 1 / 2, 1 / 2]},
              'bcc': {'G': [0, 0, 0],
                      'H': [1 / 2, -1 / 2, 1 / 2],
                      'N': [0, 0, 1 / 2],
                      'P': [1 / 4, 1 / 4, 1 / 4]},
              'hex': {'G': [0, 0, 0],
                      'M': [0, 1 / 2, 0],
                      'K': [-1 / 3, 1 / 3, 0],
                      'A': [0, 0, 1 / 2],
                      'L': [0, 1 / 2, 1 / 2],
                      'H': [-1 / 3, 1 / 3, 1 / 2]},
              'tet': {'G': [0, 0, 0],
                      'X': [1 / 2, 0, 0],
                      'M': [1 / 2, 1 / 2, 0],
                      'Z': [0, 0, 1 / 2],
                      'R': [1 / 2, 0, 1 / 2],
                      'A': [1 / 2, 1 / 2, 1 / 2]},
              'orc': {'G': [0, 0, 0],
                      'R': [1 / 2, 1 / 2, 1 / 2],
                      'S': [1 / 2, 1 / 2, 0],
                      'T': [0, 1 / 2, 1 / 2],
                      'U': [1 / 2, 0, 1 / 2],
                      'X': [1 / 2, 0, 0],
                      'Y': [0, 1 / 2, 0],
                      'Z': [0, 0, 1 / 2]}}


class SimpleBravaisLattice(BravaisLattice):
    """Special implementation for cases with only one variant."""
    special_point_names = None  # Autoinitialized by @bravais decorator
    special_points = None
    special_paths = None

    def _special_points(self, **kwargs):
        return self.special_points

    def _variant_name(self, **kwargs):
        return self.variant_names[0]


class Variant:
    variant_desc = """\
Variant name: {name}
  Special point names: {special_point_names}
  Special paths: {special_paths}
"""

    def __init__(self, name, special_point_names, special_paths):
        self.name = name
        self.special_point_names = special_point_names
        self.special_paths = special_paths

    def __str__(self):
        return self.variant_desc.format(**vars(self))


bravais_names = []
bravais_lattices = {}

def bravais(longname, parameters, variants):

    def decorate(cls):
        btype = cls.__name__
        cls.type = btype
        cls.name = longname
        cls.parameters = tuple(parameters)
        cls.variant_names = []
        cls.variants = {}

        for name, special_point_names, special_paths in variants:
            special_paths = special_paths.replace(',', ' ')
            special_paths = tuple(special_paths.split())
            cls.variant_names.append(name)
            cls.variants[name] = Variant(name, special_point_names,
                                         special_paths)

        if len(variants) == 1:
            # Only one variant.  We define the special points/paths statically:
            variant = cls.variants[cls.variant_names[0]]
            cls.special_point_names = variant.special_point_names
            cls.special_paths = variant.special_paths

            assert cls.type.isupper()
            lowername = cls.type.lower()
            from ase.dft.kpoints import ibz_points
            name2name = {'cub': 'cubic',
                         'fcc': 'fcc',
                         'bcc': 'bcc',
                         'tet': 'tetragonal',
                         'orc': 'orthorhombic',
                         'hex': 'hexagonal'}
            if lowername in name2name:
                pointinfo = ibz_points[name2name[lowername]]
                points = []
                for name in cls.special_point_names:
                    points.append(pointinfo[name])
                cls.special_points = np.array(points)

        # Register in global list and dictionary
        bravais_names.append(btype)
        bravais_lattices[btype] = cls
        return cls

    return decorate


class Cubic(SimpleBravaisLattice):
    """Abstract class for cubic lattices."""
    def __init__(self, a):
        SimpleBravaisLattice.__init__(self, a=a)

@bravais('cubic', 'a',
         [['CUB1', 'GXRM', 'GXMGRX MR']])
class CUB(Cubic):
    def _cell(self, a):
        return a * np.eye(3)

@bravais('face-centered cubic', 'a',
         [['FCC1', 'GKLUWX', 'GXWKGLUWLK UX']])
class FCC(Cubic):
    def _cell(self, a):
        return 0.5 * np.array([[0., a, a], [a, 0, a], [a, a, 0]])

@bravais('body-centered cubic', 'a',
         [['BCC1', 'GHPN', 'GHNGPH PN']])
class BCC(Cubic):
    def _cell(self, a):
        return 0.5 * np.array([[-a, a, a], [a, -a, a], [a, a, -a]])

@bravais('tetragonal', 'ac',
         [['TET1', 'GAMRXZ', 'GXMGZRAZXR MA']])
class TET(SimpleBravaisLattice):
    def __init__(self, a, c):
        SimpleBravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        return np.diag(np.array([a, a, c]))

@bravais('body-centered tetragonal', 'ac',
         [['BCT1', 'GMNPXSS1', 'GXMGSPNS1M XP'],
          ['BCT2', 'GNPSS1XYY1Z', 'GXYSGZS1NPY1Z XP']])
class BCT(BravaisLattice):
    def __init__(self, a, c):
        BravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        return 0.5 * np.array([[-a, a, c], [a, -a, c], [a, a, -c]])

    def _variant_name(self, a, c):
        return 'BCT1' if c < a else 'BCT2'

    def _special_points(self, a, c, variant):
        a2 = a * a
        c2 = c * c
        eta = .25 * (1 + c2 / a2)

        assert variant.name in self.variants

        if variant.name == 'BCT1':
            points = [[0,0,0],
                      [-.5, .5, .5],
                      [0.,.5,0.],
                      [.25, .25, .25],
                      [0.,0.,.5],
                      [eta,eta,-eta],
                      [-eta,1-eta,eta]]
        else:
            zeta = 0.5 * a2 / c2
            points = [[0.,.0,0.],
                      [0.,.5,0.],
                      [.25,.25,.25],
                      [-eta,eta,eta],
                      [eta,1-eta,-eta],
                      [0.,0.,.5],
                      [-zeta,zeta,.5],
                      [.5,.5,-zeta],
                      [.5,.5,-.5]]
        return points

class Orthorhombic(BravaisLattice):
    """Abstract class for orthorhombic types."""
    def __init__(self, a, b, c):
        BravaisLattice.__init__(self, a=a, b=b, c=c)

@bravais('orthorhombic', 'abc',
         [['ORC1', 'GRSTUXYZ', 'GXSYGZURTZ YT UX SR']])
class ORC(Orthorhombic, SimpleBravaisLattice):  # FIXME stupid diamond problem
    def _cell(self, a, b, c):
        return np.diag([a, b, c]).astype(float)

@bravais('face-centered orthorhombic', 'abc',
         [['ORCF1', 'GAA1LTXX1YZ', 'GYTZGXA1Y TX1 XAZ LG'],
          ['ORCF2', 'GCC1DD1LHH1XYZ', 'GYCDXGZD1HC C1Z XH1 HY LG'],
          ['ORCF3', 'GAA1LTXX1YZ', 'GYTZGXA1Y TX1 XAZ LG']])  # same as orcf1
class ORCF(Orthorhombic):
    def _cell(self, a, b, c):
        return 0.5 * np.array([[0, b, c], [a, 0, c], [a, b, 0]])

    def _special_points(self, a, b, c, variant):
        a2 = a * a
        b2 = b * b
        c2 = c * c
        xminus = 0.25 * (1 + a2 / b2 - a2 / c2)
        xplus = 0.25 * (1 + a2 / b2 + a2 / c2)

        variant = get_orcf_variant(a, b, c)

        if variant == 1 or variant == 3:
            zeta = xminus
            eta = xplus

            points = [[0, 0, 0],
                      [.5, .5 + zeta, zeta],
                      [.5, .5 - zeta, 1 - zeta],
                      [.5, .5, .5],
                      [1., .5, .5],
                      [0., eta, eta],
                      [1., 1 - eta, 1 - eta],
                      [.5, 0, .5],
                      [.5, .5, 0]]

        if variant == 2:
            phi = 0.25 * (1 + c2 / b2 - c2 / a2)
            delta = 0.25 * (1 + b2 / a2 - b2 / c2)
            eta = xminus

            points = [[0,0,0],
                      [.5, .5-eta, 1-eta],
                      [.5, .5+eta, eta],
                      [.5-delta, .5, 1-delta],
                      [.5+delta, .5, delta],
                      [.5, .5, .5],
                      [1-phi, .5-phi, .5],
                      [phi, .5+phi, .5],
                      [0., .5, .5],
                      [.5, 0., .5],
                      [.5, .5, 0.]]

        return points

    def _variant_name(self, a, b, c):
        # XXX
        # In general we need to remember eps, because here we need eps to
        # decide which variant we are.
        #
        # The framework must also know how to forward the eps to this function

        diff = 1.0 / (a * a) - 1.0 / (b * b) - 1.0 / (c * c)
        eps = 2e-4
        if abs(diff) < eps:
            return 'ORCF3'
        return 'ORCF1' if diff > 0 else 'ORCF2'


@bravais('body-centered orthorhombic', 'abc',
         [['ORCI1', 'GLL1L2RSTWXX1YY1Z', 'GXLTWRX1ZGYSW L1Y Y1Z']])
class ORCI(Orthorhombic):
    def _cell(self, a, b, c):
        return 0.5 * np.array([[-a, b, c], [a, -b, c], [a, b, -c]])

    def _variant_name(self, a, b, c):
        return 'ORCI1'

    def _special_points(self, a, b, c, variant):
        a2 = a**2
        b2 = b**2
        c2 = c**2

        zeta = .25 * (1 + a2 / c2)
        eta = .25 * (1 + b2 / c2)
        delta = .25 * (b2 - a2) / c2
        mu = .25 * (a2 + b2) / c2

        points = [[0.,0.,0.],
                  [-mu,-mu,.5-delta],
                  [mu, -m1, .5+delta],
                  [.5-delta, .5+delta, -mu],
                  [0,.5,0],
                  [.5,0,0],
                  [0.,0.,.5],
                  [.25,.25,.25],
                  [-zeta, zeta, zeta],
                  [zeta, 1 - zeta, -zeta],
                  [eta, -eta, eta],
                  [1 - eta, eta, -eta],
                  [.5,.5,-.5]]
        return points


@bravais('c-centered orthorhombic', 'abc',
         [['ORCC1', 'GAA1RSTXX1YZ', 'GXSRAZGYX1A1TY ZT']])
class ORCC(Orthorhombic, SimpleBravaisLattice):  # FIXME stupid diamond problem
    def _cell(self, a, b, c):
        return np.array([[0.5 * a, -0.5 * b, 0], [0.5 * a, 0.5 * b, 0],
                         [0, 0, c]])

@bravais('hexagonal', 'ac',
         [['HEX1', 'GMKALH', 'GMKGALHA LM KH']])
class HEX(SimpleBravaisLattice):
    def __init__(self, a, c):
        BravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        x = 0.5 * np.sqrt(3)
        return np.array([[0.5 * a, -x * a, 0], [0.5 * a, x * a, 0],
                         [0., 0., c]])


@bravais('rhombohedral', ('a', 'alpha'),
         [['RHL1', 'GBB21FLL1PP1P2QXZ', 'GLB1 BZGX QFP1Z LP'],
          ['RHL2', 'GFLPP1QQ1Z', 'GPZQGFP1Q1LZ']])
class RHL(BravaisLattice):
    def __init__(self, a, alpha):
         BravaisLattice.__init__(self, a=a, alpha=alpha)

    def _cell(self, a, alpha):
        alpha *= np.pi / 180
        acosa = a * np.cos(alpha)
        acosa2 = a * np.cos(0.5 * alpha)
        asina2 = a * np.sin(0.5 * alpha)
        acosfrac = acosa / acosa2
        return np.array([[acosa2, -asina2, 0], [acosa2, asina2, 0],
                         [a * acosfrac, 0, a * np.sqrt(1 - acosfrac**2)]])

    def _variant_name(self, a, alpha):
        return 'RHL1' if alpha < 90 else 'RHL2'

    def _special_points(self, a, alpha, variant):
        cosa = np.cos(alpha)
        eta = (1 + 4 * cosa) / (2 + 4 * cosa)
        nu = .75 - 0.5 * eta

        if variant.name == 'RHL1':
            points = [[0,0,0],
                      [eta,.5,1-eta],
                      [.5, 1 - eta, eta - 1],
                      [.5,.5,0],
                      [.5,0,0],
                      [0,0,-.5],
                      [eta,nu,nu],
                      [1-nu,1-nu,1-eta],
                      [nu,nu,eta-1],
                      [1-nu,nu,0],
                      [nu,0,-nu],
                      [.5,.5,.5]]
        else:
            points = [[0,0,0],
                      [.5,-.5,0],
                      [.5,0,0],
                      [1-nu,-nu,1-nu],
                      [nu,nu-1,nu-1],
                      [eta,eta,eta],
                      [1-eta,-eta,-eta],
                      [.5,-.5,.5]]
        return points

@bravais('monoclinic', ('a', 'b', 'c', 'alpha'),
         [['MCL1', 'GACDD1EHH1H2MM1M2XYY1Z', 'GYHCEM1AXH1 MDZ YD']])
class MCL(SimpleBravaisLattice):
    def __init__(self, a, b, c, alpha):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha)

    def _cell(self, a, b, c, alpha):
        alpha *= np.pi / 180
        return np.array([[a, 0, 0], [0, b, 0],
                         [0, c * np.cos(alpha), c * np.sin(alpha)]])


@bravais('c-centered monoclinic', ('a', 'b', 'c', 'alpha'),
         [['MCLC1', 'GNN1FF1F2F3II1LMXX1X2YY1Z', 'GYFLI I1ZF1 YX1 XGN MG'],
          ['MCLC2', 'GNN1FF1F2F3II1LMXX1X2YY1Z', 'GYFLI I1ZF1 NGM'],
          ['MCLC3', 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z', 'GYFHZIF1 H1Y1XGN MG'],
          ['MCLC4', 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z', 'GYFHZI H1Y1XGN MG'],
          ['MCLC5', 'GFF1F2HH1H2II1LMNN1XYY1Y2Y3Z',
           'GYFLI I1ZHF1 H1Y1XGN MG']])
class MCLC(BravaisLattice):
    def __init__(self, a, b, c, alpha):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha)

    def _cell(self, a, b, c, alpha):
        alpha *= np.pi / 180
        return np.array([[0.5 * a, 0.5 * b, 0], [-0.5 * a, 0.5 * b, 0],
                         [0, c * np.cos(alpha), c * np.sin(alpha)]])

    def _variant_name(self, a, b, c, alpha):
        #from ase.geometry.cell import mclc
        # okay, this is a bit hacky

        # We need the same parameters here as when determining the points.
        # Right now we just repeat the code:
        a2 = a * a
        b2 = b * b
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        sina2 = sina**2

        cell = self.tocell()
        lengths_angles = Cell(np.linalg.inv(cell)).cellpar()

        kgamma = lengths_angles[-1]

        eps = 2e-4  # XXX we should know precision better somehow
        # Also, we should probably not compare angles/lengths with
        # the same precision
        if abs(kgamma - 90) < eps:
            variant = 2
        elif kgamma > 90:
            variant = 1
        elif kgamma < 90:
            num = b * cosa / c + b2 * sina2 / a2
            if abs(num - 1) < eps:
                variant = 4
            elif num < 1:
                variant = 3
            else:
                variant = 5
        variant = 'MCLC' + str(variant)
        return variant

    def _special_points(self, a, b, c, alpha, variant):
        assert a <= c
        assert b <= c
        assert alpha < 90

        variant = int(variant.name[-1])

        a2 = a * a
        b2 = b * b
        # c2 = c * c
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        sina2 = sina**2

        if variant == 1 or variant == 2:
            zeta = (2 - b * cosa / c) / (4 * sina2)
            eta = 0.5 + 2 * zeta * c * cosa / b
            psi = .75 - a2 / (4 * b2 * sina * sina)
            phi = psi + (.75 - psi) * b * cosa / c

            points = [[0,0,0],
                      [.5,0,0],
                      [0,-.5,0],
                      [1-zeta,1-zeta,1-eta],
                      [zeta,zeta,eta],
                      [-zeta,-zeta,1-eta],
                      [1-zeta,-zeta,1-eta],
                      [phi,1-phi,.5],
                      [.5,.5,.5],
                      [.5,0,.5],
                      [1-psi,psi-1,0],
                      [psi,1-psi,0],
                      [.5,.5,0],
                      [-.5,-.5,0],
                      [0,0,.5]]
        elif variant == 3 or variant == 4:
            mu = .25 * (1 + b2 / a2)
            delta = b * c * cosa / (2  * a2)
            zeta = mu - 0.25 + (1 - b * cosa / c) / (4 * sina2)
            eta = 0.5 + 2 * zeta * c * cosa / b
            phi = 1 + zeta - 2 * mu
            psi = eta - 2 * delta

            points = [[0,0,0],
                      [1-phi,1-phi,1-psi],
                      [phi,phi-1,psi],
                      [1-phi,-phi,1-psi],
                      [zeta,zeta,eta],
                      [1-zeta,-zeta,1-eta],
                      [-zeta,-zeta,1-eta],
                      [.5,-.5,.5],
                      [.5,-.5,.5],
                      [.5,0,.5],
                      [.5,0,0],
                      [0,-.5,0],
                      [.5,-.5,0],
                      [mu,mu,delta],
                      [1-mu,-mu,-delta],
                      [-mu,-mu,-delta],
                      [mu,mu-1,delta],
                      [0,0,.5]]
        elif variant == 5:
            zeta = .25 * (b2 / a2 + (1 - b * cosa / c) / sina2)
            eta = 0.5 + 2 * zeta * c * cosa / b
            mu = .5 * eta + b2 / (4 * a2) - b * c * cosa / (2 * a2)
            nu = 2 * mu - zeta
            omega = (4 * nu - 1 - b2 * sina2 / a2) * c / (2 * b * cosa)
            delta = zeta * c * cosa / b + omega / 2 - .25
            rho = 1 - zeta * a2 / b2

            points = [[0,0,0],
                      [nu,nu,omega],
                      [1-nu,1-nu,1-omega],
                      [nu,nu-1,omega],
                      [zeta,zeta,eta],
                      [1-zeta,-zeta,1-eta],
                      [-zeta,-zeta,1-eta],
                      [rho,1-rho,.5],
                      [1-rho,rho-1,.5],
                      [.5,.5,.5],
                      [.5,0,.5],
                      [.5,0,0],
                      [0,-.5,0],
                      [.5,-.5,0],
                      [mu,mu,delta],
                      [1-mu,-mu,-delta],
                      [mu,mu-1,delta],
                      [0,0,.5]]

        return points


@bravais('trigonal', ('a', 'b', 'c', 'alpha', 'beta', 'gamma'),
         [['TRI1a', 'GLMNRXYZ', 'XGY LGZ NGM RG'],  # XXX labels, paths
          ['TRI2a', 'GLMNRXYZ', 'XGY LGZ NGM RG'],  # are all the same.
          ['TRI1b', 'GLMNRXYZ', 'XGY LGZ NGM RG'],
          ['TRI2b', 'GLMNRXYZ', 'XGY LGZ NGM RG']])


class TRI(BravaisLattice):
    def __init__(self, a, b, c, alpha, beta, gamma):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha, beta=beta,
                                gamma=gamma)

    def _cell(self, a, b, c, alpha, beta, gamma):
        alpha, beta, gamma = np.array([alpha, beta, gamma]) * (np.pi / 180)
        singamma = np.sin(gamma)
        cosgamma = np.cos(gamma)
        cosbeta = np.cos(beta)
        cosalpha = np.cos(alpha)
        a3x = c * cosbeta
        a3y = c / singamma * (cosalpha - cosbeta * cosgamma)
        a3z = c / singamma * np.sqrt(singamma**2 - cosalpha**2 - cosbeta**2
                                     + 2 * cosalpha * cosbeta * cosgamma)
        return np.array([[a, 0, 0], [b * cosgamma, b * singamma, 0],
                         [a3x, a3y, a3z]])

    def _variant_name(self, a, b, c, alpha, beta, gamma):
        c = Cell.new([a, b, c, alpha, beta, gamma])
        #from ase.geometry.cell import get_cell_lengths_and_angles
        #ka, kb, kc, kalpha, kbeta, kgamma = get_cell_lengths_and_angles
        #reciprocal)
        (ka, kb, kc, kalpha, kbeta,
         kgamma) = Cell(c.reciprocal()).cellpar()
        # lengths = np.array([ka, kb, kc])
        angles = np.array([kalpha, kbeta, kgamma])

        eps = 2e-4  # XXX must support variable eps
        if abs(kgamma - 90) < eps:
            if kalpha > 90 and kbeta > 90:
                var = '2a'
            elif kalpha < 90 and kbeta < 90:
                var = '2b'
            else:
                bad
        elif all(angles > 90) and kgamma < min(kalpha, kbeta):
            var = '1a'
        elif all(angles < 90) and kgamma > max(kalpha, kbeta):
            var = '1b'
        return 'TRI' + var

    def _special_points(self, a, b, c, alpha, beta, gamma, variant):
        # (None of the points actually depend on any parameters)
        # (We should store the points openly on the variant objects)
        if var == 'TRI1a' or var == 'TRI2a':
            points = [[0.,0.,0.],
                      [.5,.5,0],
                      [0,.5,.5],
                      [.5,0,.5],
                      [.5,.5,.5],
                      [.5,0,0],
                      [0,.5,0],
                      [0,0,.5]]
        else:
            points = [[0,0,0],
                      [.5,-.5,0],
                      [0,0,.5],
                      [-.5,-.5,.5],
                      [0,-.5,.5],
                      [0,-0.5,0],
                      [.5,0,0],
                      [-.5,0,.5]]

        return points


def crystal_structure_from_cell(cell, eps=2e-4, niggli_reduce=True):
    """Return the crystal structure as a string calculated from the cell.

    Supply a cell (from atoms.get_cell()) and get a string representing
    the crystal structure returned. Works exactly the opposite
    way as ase.dft.kpoints.get_special_points().

    Parameters:

    cell : numpy.array or list
        An array like atoms.get_cell()

    Returns:

    crystal structure : str
        'cubic', 'fcc', 'bcc', 'tetragonal', 'orthorhombic',
        'hexagonal' or 'monoclinic'
    """
    cellpar = cell_to_cellpar(cell)
    abc = cellpar[:3]
    angles = cellpar[3:] / 180 * pi
    a, b, c = abc
    alpha, beta, gamma = angles

    if abc.ptp() < eps and abs(angles - pi / 2).max() < eps:
        return 'cubic'
    elif abc.ptp() < eps and abs(angles - pi / 3).max() < eps:
        return 'fcc'
    elif abc.ptp() < eps and abs(angles - np.arccos(-1 / 3)).max() < eps:
        return 'bcc'
    elif abs(a - b) < eps and abs(angles - pi / 2).max() < eps:
        return 'tetragonal'
    elif abs(angles - pi / 2).max() < eps:
        return 'orthorhombic'
    elif (abs(a - b) < eps and
          (abs(gamma - pi / 3 * 2) < eps or abs(gamma - pi / 3) < eps) and
          abs(angles[:2] - pi / 2).max() < eps):
        return 'hexagonal'
    elif (abs(angles - pi / 2) > eps).sum() == 1:
        return 'monoclinic'
    elif (abc.ptp() < eps and angles.ptp() < eps and
          np.abs(angles).max() < pi / 2):
        return 'rhombohedral type 1'
    elif (abc.ptp() < eps and angles.ptp() < eps and
          np.abs(angles).max() > pi / 2):
        return 'rhombohedral type 2'
    else:
        if niggli_reduce:
            from ase.build.tools import niggli_reduce_cell
            cell, _ = niggli_reduce_cell(cell)
            return crystal_structure_from_cell(cell, niggli_reduce=False)
        raise ValueError('Cannot find crystal structure')


def complete_cell(cell):
    """Calculate complete cell with missing lattice vectors.

    Returns a new 3x3 ndarray.
    """

    cell = np.array(cell, dtype=float)
    missing = np.nonzero(~cell.any(axis=1))[0]

    if len(missing) == 3:
        cell.flat[::4] = 1.0
    if len(missing) == 2:
        # Must decide two vectors:
        i = 3 - missing.sum()
        assert abs(cell[i, missing]).max() < 1e-16, "Don't do that"
        cell[missing, missing] = 1.0
    elif len(missing) == 1:
        i = missing[0]
        cell[i] = np.cross(cell[i - 2], cell[i - 1])
        cell[i] /= np.linalg.norm(cell[i])

    return cell


def is_orthorhombic(cell):
    """Check that cell only has stuff in the diagonal."""
    return not (np.flatnonzero(cell) % 4).any()


def orthorhombic(cell):
    """Return cell as three box dimensions or raise ValueError."""
    if not is_orthorhombic(cell):
        raise ValueError('Not orthorhombic')
    return cell.diagonal().copy()



def get_bravais_lattice1(uc, eps=2e-4):
    if np.linalg.det(uc.array) < 0:
        raise ValueError('Cell should be right-handed')

    cellpar = uc.cellpar()
    ABC = cellpar[:3]
    angles = cellpar[3:]
    A, B, C, alpha, beta, gamma = cellpar

    def categorize_differences(numbers):
        a, b, c = numbers
        eq = [abs(b - c) < eps, abs(c - a) < eps, abs(a - b) < eps]
        neq = sum(eq)

        all_equal = neq == 3
        all_different = neq == 0
        funny_direction = np.argmax(eq) if neq == 1 else None
        assert neq != 2
        return all_equal, all_different, funny_direction

    (all_lengths_equal, all_lengths_different,
     unequal_length_dir) = categorize_differences(ABC)

    (all_angles_equal, all_angles_different,
     unequal_angle_dir) = categorize_differences(angles)

    def check(f, *args, **kwargs):
        axis = kwargs.pop('axis', 0)
        cell = f(*args, **kwargs).tocell()
        mycellpar = Cell(cell).cellpar()
        permutation = (np.arange(-3, 0) + axis) % 3
        mycellpar = mycellpar.reshape(2, 3)[:, permutation].ravel()
        if np.allclose(mycellpar, cellpar):
            # Return bravais function as well as the bravais parameters
            # that would reproduce the cell
            d = dict(zip(f.parameters, args))
            d.update(kwargs)
            #if axis:
            #    d['cycle'] = axis
            return f, d

    _c = uc.array
    BC_CA_AB = np.array([np.vdot(_c[1], _c[2]),
                         np.vdot(_c[2], _c[0]),
                         np.vdot(_c[0], _c[1])])

    _, _, unequal_scalarprod_dir = categorize_differences(BC_CA_AB)

    def allclose(a, b):
        return np.allclose(a, b, atol=eps)

    if all_lengths_equal:
        if allclose(angles, 90):
            return check(CUB, A)
        if allclose(angles, 60):
            return check(FCC, np.sqrt(2) * A)
        if allclose(angles, np.arccos(-1 / 3) * 180 / np.pi):
            return check(BCC, 2.0 * A / np.sqrt(3))

    if all_lengths_equal and unequal_angle_dir is not None:
        x = BC_CA_AB[unequal_angle_dir]
        y = BC_CA_AB[(unequal_angle_dir + 1) % 3]

        if x < 0:
            c = 2.0 * np.sqrt(-y)
            a = np.sqrt(2.0 * A**2 - 0.5 * c**2)
            obj = check(BCT, a, c, axis=-unequal_angle_dir + 2)
            if obj:
                return obj

    if (unequal_angle_dir is not None
          and abs(angles[unequal_angle_dir] - 120) < eps
          and abs(angles[unequal_angle_dir - 1] - 90) < eps):
        a2 = -2 * BC_CA_AB[unequal_scalarprod_dir]
        c = ABC[unequal_scalarprod_dir]
        assert a2 > 0
        return check(HEX, np.sqrt(a2), c, axis=-unequal_scalarprod_dir + 2)

    if allclose(angles, 90) and unequal_length_dir is not None:
        a = ABC[unequal_length_dir - 1]
        c = ABC[unequal_length_dir]
        return check(TET, a, c, axis=-unequal_length_dir + 2)

    if unequal_length_dir is not None:
        X = ABC[unequal_length_dir - 1]**2
        Y = BC_CA_AB[unequal_length_dir]
        c = ABC[unequal_length_dir]
        a = np.sqrt(2 * (X + Y))
        b = np.sqrt(2 * (X - Y))
        obj = check(ORCC, a, b, c, axis=2 - unequal_length_dir)
        if obj:
            return obj

    if allclose(angles, 90) and all_lengths_different:
        return check(ORC, A, B, C)

    if all_lengths_different:
        obj = check(ORCF, *(2 * np.sqrt(BC_CA_AB)))
        if obj:
            return obj

    if all_lengths_equal:
        dims2 = -2 * np.array([BC_CA_AB[1] + BC_CA_AB[2],
                               BC_CA_AB[2] + BC_CA_AB[0],
                               BC_CA_AB[0] + BC_CA_AB[1]])
        if all(dims2 > 0):
            dims = np.sqrt(dims2)
            obj = check(ORCI, *dims)
            if obj:
                return obj

    if all_lengths_equal:
        cosa = BC_CA_AB[0] / A**2
        alpha = np.arccos(cosa) * 180 / np.pi
        obj = check(RHL, A, alpha)
        if obj:
            return obj

    if all_lengths_different and unequal_scalarprod_dir is not None:
        alpha = angles[unequal_scalarprod_dir]
        abc = ABC[np.arange(-3, 0) + unequal_scalarprod_dir]
        obj = check(MCL, *abc, alpha=alpha, axis=-unequal_scalarprod_dir)
        if obj:
            return obj

    if unequal_length_dir is not None:
        c = ABC[unequal_length_dir]
        L = ABC[unequal_length_dir - 1]
        b = np.sqrt(2 * (L**2 + BC_CA_AB[unequal_length_dir]))
        a = np.sqrt(4 * L**2 - b**2)
        cosa = 2 * BC_CA_AB[unequal_length_dir - 1] / (b * c)
        alpha = np.arccos(cosa) * 180 / np.pi
        obj = check(MCLC, a, b, c, alpha, axis=-unequal_length_dir + 2)
        if obj:
            return obj

    obj = check(TRI, A, B, C, *angles)
    if obj:
        # Should always be true
        return obj

    raise RuntimeError('Cannot recognize cell at all somehow!')

def get_bravais_lattice(uc, eps=2e-4):
    bravaisclass, parameters = get_bravais_lattice1(uc, eps=eps)
    # XXX MCL does not get alpha.  Only a, b, c
    bravais = bravaisclass(**parameters)
    return bravais
