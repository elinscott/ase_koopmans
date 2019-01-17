from __future__ import division
from abc import abstractmethod, ABC
import re

import numpy as np

from ase.geometry.cell import Cell

_degrees = np.pi / 180


class BandPath:
    def __init__(self,  icell, coords, xvalues, labels=None, special_coords=None):
        if labels is None and special_coords is None:
            labels = ''
            special_coords = np.empty((0, 3))
        assert len(xvalues) == len(coords)
        self.xvalues = xvalues
        assert len(labels) == len(special_coords)  # Maybe include ','?
        self.labels = labels
        assert special_coords.shape == (labels, 3)
        self.special_coords = special_coords

        self.icell = icell
        self.coords = coords

    def _scale(self, coords):
        return np.dot(coords, self.icell)

    @property
    def array(self):
        return self._scale(self.coords)


class BravaisLattice(ABC):
    # These parameters can be set by the @bravais decorator for a subclass.
    # (We could also use metaclasses to do this, but that's more abstract)
    type = None  # e.g. 'CUB', 'BCT', 'ORCF', ...
    name = None  # e.g. 'cubic', 'body-centered tetragonal', ...
    parameters = None  # e.g. ('a', 'c')
    variants = None  # e.g. {'BCT1': <variant object>,
    #                        'BCT2': <variant object>}

    def __init__(self, **kwargs):
        p = {}
        eps = kwargs.pop('eps', 2e-4)
        for k, v in kwargs.items():
            p[k] = float(v)
        self._parameters = p
        self._eps = eps
        self._variant = self.get_variant()

    @property
    def variant(self):
        return self._variant

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

    @property
    def special_path(self):
        return self.variant.special_path

    def get_special_points(self):
        variant = self._variant
        points = self._special_points(variant=variant, **self._parameters)
        assert len(points) == len(self.kpoint_labels)
        return np.array(points)

    # XXX which should be the standard way of doing this?
    # Array or dict?
    def get_special_point_dict(self):
        labels = self.kpoint_labels
        points = self.get_special_points()
        return dict(zip(labels, points))

    def get_variant(self):
        name = self._variant_name(**self._parameters)
        return self.variants[name]

    @property
    def kpoint_labels(self):
        labels = re.findall(r'[A-Z]\d?', self._variant.special_point_names)
        return labels

    def plot_bz(self, path=None, **plotkwargs):
        from ase.dft.bz import bz3d_plot

        coords = self.get_special_point_dict()

        if path is None:
            isolated_points = ','.join(coords)
            path = self.special_path + ',' + isolated_points
            # (Isolated points are normally plotted twice since they are
            #  also part of the special path.)

        cell = self.tocell()
        icell = cell.reciprocal()

        paths = []

        if path:
            for path0 in path.split(','):
                if not re.match(r'([A-Z]\d?)+', path0):
                    raise ValueError('Invalid path string: {}'
                                     .format(repr(path0)))
                path0 = re.findall(r'[A-Z]\d?', path0)
                thecoords = [coords[label] for label in path0]
                abscoords = np.dot(thecoords, icell)
                paths.append((path0, abscoords))
        else:
            paths = None

        kw = {'vectors': True}
        kw.update(plotkwargs)

        return bz3d_plot(cell, paths=paths, **kw)

    def bandpath(self, path=None, npoints=50):
        # npoints should depend on the length of the path
        if path is None:
            path = self.variant.special_path

        pathcoords = self._resolve_kpt_path_string(path)

        from ase.dft.kpoints import paths2kpts
        cell = self.tocell()
        icell = cell.reciprocal()
        kpts, x, X = paths2kpts(pathcoords, icell, npoints)

        return BandPath(icell, kpts, x, labels=path, special_coords=X)

    def _resolve_kpt_path_string(self, path):
        from ase.dft.kpoints import parse_path_string
        paths = parse_path_string(path)
        special_point_coords = self.get_special_point_dict()
        return [np.array([special_point_coords[sym] for sym in subpath])
                for subpath in paths]

    @abstractmethod
    def _cell(self, **kwargs):
        """Return a Cell object from this Bravais lattice.

        Arguments are the dictionary of Bravais parameters."""
        pass

    @abstractmethod
    def _special_points(self, **kwargs):
        """Return the special point coordinates as an npoints x 3 sequence.

        Subclasses typically return a nested list.

        Ordering must be same as kpoint labels.

        Arguments are the dictionary of Bravais parameters and the variant."""
        pass

    @abstractmethod
    def _variant_name(self, **kwargs):
        """Return the name (e.g. ORCF3) of variant.

        Arguments will be the dictionary of Bravais parameters."""
        pass

    def __repr__(self):
        par = ', '.join('{}={}'.format(k, v)
                        for k, v in self._parameters.items())
        return '{}({})'.format(self.type, par)

    def __str__(self):
        points = self.get_special_points()
        labels = self.kpoint_labels

        coordstring = '\n'.join(['    {:2s} {:7.4f} {:7.4f} {:7.4f}'
                                 .format(label, *point)
                                 for label, point
                                 in zip(labels, points)])

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
        desc = """\
Lattice name: {type}
  Long name: {name}
  Parameters: {parameters}
""".format(**vars(cls))

        chunks = [desc]
        for name in cls.variant_names:
            var = cls.variants[name]
            txt = str(var)
            lines = ['  ' + L for L in txt.splitlines()]
            lines.append('')
            chunks.extend(lines)

        return '\n'.join(chunks)


class SimpleBravaisLattice(BravaisLattice):
    """Special implementation for cases with only one variant."""
    special_point_names = None  # These are initialized by @bravais decorator
    special_points = None
    special_path = None

    def _special_points(self, **kwargs):
        assert self.special_points is not None, self.__class__
        return self.special_points

    def _variant_name(self, **kwargs):
        assert len(self.variants) == 1
        return self.variant_names[0]


class Variant:
    variant_desc = """\
Variant name: {name}
  Special point names: {special_point_names}
  Default path: {special_path}
"""

    def __init__(self, name, special_point_names, special_path):
        self.name = name
        self.special_point_names = special_point_names
        self.special_path = special_path

    def __str__(self):
        return self.variant_desc.format(**vars(self))


bravais_names = []
bravais_lattices = {}

def bravais(longname, parameters, variants):
    """Decorator for Bravais lattice classes.

    This sets a number of class variables and processes the information
    about different variants into a list of Variant objects."""

    def decorate(cls):
        btype = cls.__name__
        cls.type = btype
        cls.name = longname
        cls.parameters = tuple(parameters)
        cls.variant_names = []
        cls.variants = {}

        for name, special_point_names, special_path in variants:
            cls.variant_names.append(name)
            cls.variants[name] = Variant(name, special_point_names,
                                         special_path)

        if len(variants) == 1:
            # Only one variant.  We define the special points/paths statically:
            variant = cls.variants[cls.variant_names[0]]
            cls.special_point_names = variant.special_point_names
            cls.special_path = variant.special_path

            assert cls.type.isupper()
            lowername = cls.type.lower()
            from ase.dft.kpoints import ibz_points
            name2name = {'cub': 'cubic',
                         'fcc': 'fcc',
                         'bcc': 'bcc',
                         'tet': 'tetragonal',
                         'orc': 'orthorhombic'}
            pointinfo = None
            if lowername in name2name:
                pointinfo = ibz_points[name2name[lowername]]
            elif lowername == 'hex':
                pointinfo = {'Gamma': [0, 0, 0],
                              'M': [0, 1 / 2, 0],
                              'K': [1 / 3, 1 / 3, 0],  # !
                              'A': [0, 0, 1 / 2],
                              'L': [0, 1 / 2, 1 / 2],
                              'H': [1 / 3, 1 / 3, 1 / 2]} # !

            if pointinfo is not None:
                points = []
                for name in cls.special_point_names:
                    if name == 'G':
                        name = 'Gamma'
                    points.append(pointinfo[name])
                cls.special_points = np.array(points)

        # Register in global list and dictionary
        bravais_names.append(btype)
        bravais_lattices[btype] = cls
        return cls

    return decorate


class UnconventionalLattice(ValueError):
    pass


class Cubic(SimpleBravaisLattice):
    """Abstract class for cubic lattices."""
    def __init__(self, a):
        SimpleBravaisLattice.__init__(self, a=a)

@bravais('cubic', 'a',
         [['CUB1', 'GXRM', 'GXMGRX,MR']])
class CUB(Cubic):
    def _cell(self, a):
        return a * np.eye(3)

@bravais('face-centered cubic', 'a',
         [['FCC1', 'GKLUWX', 'GXWKGLUWLK,UX']])
class FCC(Cubic):
    def _cell(self, a):
        return 0.5 * np.array([[0., a, a], [a, 0, a], [a, a, 0]])

@bravais('body-centered cubic', 'a',
         [['BCC1', 'GHPN', 'GHNGPH,PN']])
class BCC(Cubic):
    def _cell(self, a):
        return 0.5 * np.array([[-a, a, a], [a, -a, a], [a, a, -a]])

@bravais('tetragonal', 'ac',
         [['TET1', 'GAMRXZ', 'GXMGZRAZ,XR,MA']])
class TET(SimpleBravaisLattice):
    def __init__(self, a, c):
        SimpleBravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        return np.diag(np.array([a, a, c]))

# XXX in BCT2 we use S for Sigma.
# Also in other places I think
@bravais('body-centered tetragonal', 'ac',
         [['BCT1', 'GMNPXZZ1', 'GXMGZPNZ1M,XP'],
          ['BCT2', 'GNPSS1XYY1Z', 'GXYSGZS1NPY1Z,XP']])
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

        assert variant.name in self.variants

        if variant.name == 'BCT1':
            eta = .25 * (1 + c2 / a2)
            points = [[0,0,0],
                      [-.5, .5, .5],
                      [0.,.5,0.],
                      [.25, .25, .25],
                      [0.,0.,.5],
                      [eta,eta,-eta],
                      [-eta,1-eta,eta]]
        else:
            eta = .25 * (1 + a2 / c2)  # Not same eta as BCT1!
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


def check_orc(a, b, c):
    if not a < b < c:
        raise UnconventionalLattice('Expected a < b < c, got {}, {}, {}'
                                    .format(a, b, c))


class Orthorhombic(BravaisLattice):
    """Abstract class for orthorhombic types."""
    def __init__(self, a, b, c):
        BravaisLattice.__init__(self, a=a, b=b, c=c)


@bravais('orthorhombic', 'abc',
         [['ORC1', 'GRSTUXYZ', 'GXSYGZURTZ,YT,UX,SR']])
class ORC(Orthorhombic, SimpleBravaisLattice):  # FIXME stupid diamond problem
    def _cell(self, a, b, c):
        return np.diag([a, b, c]).astype(float)


@bravais('face-centered orthorhombic', 'abc',
         [['ORCF1', 'GAA1LTXX1YZ', 'GYTZGXA1Y,TX1,XAZ,LG'],
          ['ORCF2', 'GCC1DD1LHH1XYZ', 'GYCDXGZD1HC,C1Z,XH1,HY,LG'],
          ['ORCF3', 'GAA1LTXX1YZ', 'GYTZGXA1Y,XAZ,LG']])
class ORCF(Orthorhombic):
    def _cell(self, a, b, c):
        return 0.5 * np.array([[0, b, c], [a, 0, c], [a, b, 0]])

    def _special_points(self, a, b, c, variant):
        a2 = a * a
        b2 = b * b
        c2 = c * c
        xminus = 0.25 * (1 + a2 / b2 - a2 / c2)
        xplus = 0.25 * (1 + a2 / b2 + a2 / c2)

        if variant.name == 'ORCF1' or variant.name == 'ORCF3':
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
        else:
            assert variant.name == 'ORCF2'
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
        check_orc(a, b, c)

        diff = 1.0 / (a * a) - 1.0 / (b * b) - 1.0 / (c * c)
        if abs(diff) < self._eps:
            return 'ORCF3'
        return 'ORCF1' if diff > 0 else 'ORCF2'

@bravais('body-centered orthorhombic', 'abc',
         [['ORCI1', 'GLL1L2RSTWXX1YY1Z', 'GXLTWRX1ZGYSW,L1Y,Y1Z']])
class ORCI(Orthorhombic):
    def _cell(self, a, b, c):
        return 0.5 * np.array([[-a, b, c], [a, -b, c], [a, b, -c]])

    def _variant_name(self, a, b, c):
        check_orc(a, b, c)
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
                  [-mu,mu,.5-delta],
                  [mu, -mu, .5+delta],
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
         [['ORCC1', 'GAA1RSTXX1YZ', 'GXSRAZGYX1A1TY,ZT']])
class ORCC(Orthorhombic, SimpleBravaisLattice):  # FIXME stupid diamond problem
    def _cell(self, a, b, c):
        return np.array([[0.5 * a, -0.5 * b, 0], [0.5 * a, 0.5 * b, 0],
                         [0, 0, c]])

    def _special_points(self, a, b, c, variant):
        zeta = .25 * (1 + a * a / (b * b))
        points = [[0,0,0],
                  [zeta,zeta,.5],
                  [-zeta,1-zeta,.5],
                  [0,.5,.5],
                  [0,.5,0],
                  [-.5,.5,.5],
                  [zeta,zeta,0],
                  [-zeta,1-zeta,0],
                  [-.5,.5,0],
                  [0,0,.5]]
        return points


    def _variant_name(self, a, b, c):
        if not a < b:
            raise UnconventionalLattice('Expected a < b, got a={}, b={}'
                                        .format(a, b))
        return SimpleBravaisLattice._variant_name(self)

@bravais('hexagonal', 'ac',
         [['HEX1', 'GMKALH', 'GMKGALHA,LM,KH']])
class HEX(SimpleBravaisLattice):
    def __init__(self, a, c):
        BravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        x = 0.5 * np.sqrt(3)
        return np.array([[0.5 * a, -x * a, 0], [0.5 * a, x * a, 0],
                         [0., 0., c]])


@bravais('rhombohedral', ('a', 'alpha'),
         [['RHL1', 'GBB1FLL1PP1P2QXZ', 'GLB1,BZGX,QFP1Z,LP'],
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
                         [a * acosfrac, 0, a * (1 - acosfrac**2)**.5]])

    def _variant_name(self, a, alpha):
        return 'RHL1' if alpha < 90 else 'RHL2'

    def _special_points(self, a, alpha, variant):
        if variant.name == 'RHL1':
            cosa = np.cos(alpha * _degrees)
            eta = (1 + 4 * cosa) / (2 + 4 * cosa)
            nu = .75 - 0.5 * eta
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
            eta = 1 / (2 * np.tan(alpha * _degrees / 2)**2)
            nu = .75 - 0.5 * eta
            points = [[0,0,0],
                      [.5,-.5,0],
                      [.5,0,0],
                      [1-nu,-nu,1-nu],
                      [nu,nu-1,nu-1],
                      [eta,eta,eta],
                      [1-eta,-eta,-eta],
                      [.5,-.5,.5]]
        return points


def check_mcl(a, b, c, alpha):
    if not (a <= c and b <= c and alpha < 90):
        raise UnconventionalLattice('Expected a <= c, b <= c, alpha < 90')


@bravais('monoclinic', ('a', 'b', 'c', 'alpha'),
         [['MCL1', 'GACDD1EHH1H2MM1M2XYY1Z', 'GYHCEM1AXH1,MDZ,YD']])
class MCL(BravaisLattice):
    def __init__(self, a, b, c, alpha):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha)

    def _cell(self, a, b, c, alpha):
        alpha *= np.pi / 180
        return np.array([[a, 0, 0], [0, b, 0],
                         [0, c * np.cos(alpha), c * np.sin(alpha)]])

    def _special_points(self, a, b, c, alpha, variant):
        cosa = np.cos(alpha * _degrees)
        eta = (1 - b * cosa / c) / (2 * np.sin(alpha)**2)
        nu = .5 - eta * c * cosa / b

        points = [[0,0,0],
                  [.5,.5,0],
                  [0,.5,.5],
                  [.5,0,.5],
                  [.5,0,-.5],
                  [.5,.5,.5],
                  [0,eta,1-nu],
                  [0,1-eta,nu],
                  [0,eta,-nu],
                  [.5,eta,1-nu],
                  [.5,1-eta,nu],
                  [.5,eta,-nu],
                  [0,.5,0],
                  [0,0,.5],
                  [0,0,-.5],
                  [.5,0,0]]
        return points

    def _variant_name(self, a, b, c, alpha):
        check_mcl(a, b, c, alpha)
        return 'MCL1'


@bravais('c-centered monoclinic', ('a', 'b', 'c', 'alpha'),
         [['MCLC1', 'GNN1FF1F2F3II1LMXX1X2YY1Z', 'GYFLI,I1ZF1,YX1,XGN,MG'],
          ['MCLC2', 'GNN1FF1F2F3II1LMXX1X2YY1Z', 'GYFLI,I1ZF1,NGM'],
          ['MCLC3', 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z', 'GYFHZIF1,H1Y1XGN,MG'],
          ['MCLC4', 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z', 'GYFHZI,H1Y1XGN,MG'],
          ['MCLC5', 'GFF1F2HH1H2II1LMNN1XYY1Y2Y3Z',
           'GYFLI,I1ZHF1,H1Y1XGN,MG']])
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
        check_mcl(a, b, c, alpha)

        a2 = a * a
        b2 = b * b
        cosa = np.cos(alpha * _degrees)
        sina = np.sin(alpha * _degrees)
        sina2 = sina**2

        cell = self.tocell()
        lengths_angles = Cell(cell.reciprocal()).cellpar()

        kgamma = lengths_angles[-1]

        eps = self._eps
        # We should not compare angles in degrees versus lengths with
        # the same precision.
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
        variant = int(variant.name[-1])

        a2 = a * a
        b2 = b * b
        # c2 = c * c
        cosa = np.cos(alpha * _degrees)
        sina = np.sin(alpha * _degrees)
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
                      [1-phi,phi-1,.5],
                      [.5,.5,.5],
                      [.5,0,.5],
                      [1-psi,psi-1,0],
                      [psi,1-psi,0],
                      [psi-1,-psi,0],
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
                      [-mu,-mu,-delta],
                      [mu,mu-1,delta],
                      [0,0,.5]]

        return points


@bravais('trigonal', ('a', 'b', 'c', 'alpha', 'beta', 'gamma'),
         [['TRI1a', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG'],  # XXX labels, paths
          ['TRI2a', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG'],  # are all the same.
          ['TRI1b', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG'],
          ['TRI2b', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG']])
class TRI(BravaisLattice):
    def __init__(self, a, b, c, alpha, beta, gamma):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha, beta=beta,
                                gamma=gamma)

    def _cell(self, a, b, c, alpha, beta, gamma):
        alpha, beta, gamma = np.array([alpha, beta, gamma]) * (np.pi / 180)
        singamma = np.sin(gamma * _degrees)
        cosgamma = np.cos(gamma * _degrees)
        cosbeta = np.cos(beta * _degrees)
        cosalpha = np.cos(alpha * _degrees)
        a3x = c * cosbeta
        a3y = c / singamma * (cosalpha - cosbeta * cosgamma)
        a3z = c / singamma * np.sqrt(singamma**2 - cosalpha**2 - cosbeta**2
                                     + 2 * cosalpha * cosbeta * cosgamma)
        return np.array([[a, 0, 0], [b * cosgamma, b * singamma, 0],
                         [a3x, a3y, a3z]])

    def _variant_name(self, a, b, c, alpha, beta, gamma):
        c = Cell.new([a, b, c, alpha, beta, gamma])
        icellpar = Cell(c.reciprocal()).cellpar()
        kangles = kalpha, kbeta, kgamma = icellpar[3:]

        eps = self._eps
        if abs(kgamma - 90) < eps:
            if kalpha > 90 and kbeta > 90:
                var = '2a'
            elif kalpha < 90 and kbeta < 90:
                var = '2b'
            else:
                # Is this possible?  Maybe due to epsilon
                assert 0, 'unexpected combination of angles'
        elif all(kangles > 90):# and kgamma < min(kalpha, kbeta):
            var = '1a'
        elif all(kangles < 90):# and kgamma > max(kalpha, kbeta):
            var = '1b'
        else:
            raise UnconventionalLattice(
                'Reciprocal lattice has unexpected angles: kalpha={}, '
                'kbeta={}, kgamma={}'.format(kalpha, kbeta, kgamma))
        return 'TRI' + var

    def _special_points(self, a, b, c, alpha, beta, gamma, variant):
        # (None of the points actually depend on any parameters)
        # (We should store the points openly on the variant objects)
        if variant.name == 'TRI1a' or variant.name == 'TRI2a':
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


def get_bravais_lattice1(uc, eps=2e-4):
    uc2 = uc.niggli_reduce()
    if np.abs(uc2 - uc).max() > eps:
        raise ValueError('Can only get recognize Bravais lattice of '
                         'Niggli-reduced cell.')
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
        try:
            cell = f(*args, **kwargs).tocell()
        except UnconventionalLattice:
            return None
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


def _test_all_variants():
    """For testing; yield every variant of every Bravais lattice."""
    a, b, c = 3., 4., 5.
    alpha, beta, gamma = 80.0, 75.0, 65.0
    yield CUB(a)
    yield FCC(a)
    yield BCC(a)
    yield TET(a, c)
    bct1 = BCT(2 * a, c)
    bct2 = BCT(a, c)
    assert bct1.variant.name == 'BCT1'
    assert bct2.variant.name == 'BCT2'

    yield bct1
    yield bct2

    yield ORC(a, b, c)

    a0 = np.sqrt(1.0 / (1 / b**2 + 1 / c**2))
    orcf1 = ORCF(0.5 * a0, b, c)
    orcf2 = ORCF(1.2 * a0, b, c)
    orcf3 = ORCF(a0, b, c)
    assert orcf1.variant.name == 'ORCF1'
    assert orcf2.variant.name == 'ORCF2'
    assert orcf3.variant.name == 'ORCF3'
    yield orcf1
    yield orcf2
    yield orcf3

    yield ORCI(a, b, c)
    yield ORCC(a, b, c)

    yield HEX(a, c)

    rhl1 = RHL(a, alpha)
    rhl2 = RHL(a, alpha + 2.0 * (90 - alpha))
    assert rhl1.variant.name == 'RHL1'
    assert rhl2.variant.name == 'RHL2'
    yield rhl1
    yield rhl2

    yield MCL(a, b, c, alpha)

    mclc1 = MCLC(a, b, c, 80)
    mclc3 = MCLC(1.8 * a, b, c * 2, 80)
    assert mclc1.variant.name == 'MCLC1'
    assert mclc3.variant.name == 'MCLC3'

    yield mclc1
    # mclc2 has same special points as mclc1
    yield mclc3
    # mclc4 has same special points as mclc3

    mclc5 = MCLC(b, b, b, 50)
    yield mclc5

    # wtf, these are weird
    #tri1a = TRI(a*1.2,a*1.05, a*1.1, 70,80,88)
    #yield tri1a

    # XXX TODO:
    # tri1a
    # tri1b
    # tri2a
    # tri2b

def get_bravais_lattice(uc, eps=2e-4):
    bravaisclass, parameters = get_bravais_lattice1(uc, eps=eps)
    # XXX MCL does not get alpha.  Only a, b, c
    bravais = bravaisclass(**parameters)
    return bravais
