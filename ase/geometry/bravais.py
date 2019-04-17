from __future__ import division
from abc import abstractmethod, ABC

import numpy as np

from ase.geometry.cell import Cell
from ase.dft.kpoints import parse_path_string, ibz_points, BandPath


_degrees = np.pi / 180


class BravaisLattice(ABC):
    # These parameters can be set by the @bravais decorator for a subclass.
    # (We could also use metaclasses to do this, but that's more abstract)
    name = None  # e.g. 'CUB', 'BCT', 'ORCF', ...
    longname = None  # e.g. 'cubic', 'body-centered tetragonal', ...
    parameters = None  # e.g. ('a', 'c')
    variants = None  # e.g. {'BCT1': <variant object>,
    #                        'BCT2': <variant object>}

    def __init__(self, **kwargs):
        p = {}
        eps = kwargs.pop('eps', 2e-4)
        for k, v in kwargs.items():
            p[k] = float(v)
        assert set(p) == set(self.parameters)
        self._parameters = p
        self._eps = eps

        if len(self.variants) == 1:
            # If there's only one it has the same name as the lattice type
            self._variant = self.variants[self.name]
        else:
            name = self._variant_name(**self._parameters)
            self._variant = self.variants[name]

    @property
    def variant(self):
        return self._variant

    def __getattr__(self, name):
        if name in self._parameters:
            return self._parameters[name]
        return self.__getattribute__(name)  # Raises error

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

    def get_special_points_array(self):
        if self.variant.special_points is not None:
            # Fixed dictionary of special points
            d = self.get_special_points()
            labels = self.special_point_names
            assert len(d) == len(labels)
            points = np.empty((len(d), 3))
            for i, label in enumerate(labels):
                points[i] = d[label]
            return points

        # Special points depend on lattice parameters:
        points = self._special_points(variant=self.variant,
                                      **self._parameters)
        assert len(points) == len(self.special_point_names)
        return np.array(points)

    def get_special_points(self):
        if self.variant.special_points is not None:
            return self.variant.special_points

        labels = self.special_point_names
        points = self.get_special_points_array()
        return dict(zip(labels, points))

    @property
    def special_point_names(self):
        labels = parse_path_string(self._variant.special_point_names)
        assert len(labels) == 1  # list of lists
        return labels[0]

    def plot_bz(self, path=None, special_points=None, **plotkwargs):
        # Create a generic bandpath (no actual kpoints):
        bandpath = self.bandpath(path=path, special_points=special_points,
                                 npoints=0)
        return bandpath.plot(**plotkwargs)

    def bandpath(self, path=None, npoints=None, special_points=None, density=None):
        if special_points is None:
            special_points = self.get_special_points()

        if path is None:
            path = self.variant.special_path

        bandpath = BandPath(cell=self.tocell(), labelseq=path,
                            special_points=special_points)
        return bandpath.interpolate(npoints=npoints, density=density)

    @abstractmethod
    def _cell(self, **kwargs):
        """Return a Cell object from this Bravais lattice.

        Arguments are the dictionary of Bravais parameters."""
        pass

    def _special_points(self, **kwargs):
        """Return the special point coordinates as an npoints x 3 sequence.

        Subclasses typically return a nested list.

        Ordering must be same as kpoint labels.

        Arguments are the dictionary of Bravais parameters and the variant."""
        raise NotImplementedError

    def _variant_name(self, **kwargs):
        """Return the name (e.g. ORCF3) of variant.

        Arguments will be the dictionary of Bravais parameters."""
        raise NotImplementedError

    def __repr__(self):
        par = ', '.join('{}={}'.format(k, v)
                        for k, v in self._parameters.items())
        return '{}({})'.format(self.name, par)

    def __str__(self):
        points = self.get_special_points()
        labels = self.special_point_names

        coordstring = '\n'.join(['    {:2s} {:7.4f} {:7.4f} {:7.4f}'
                                 .format(label, *points[label])
                                 for label in labels])

        string = """\
{repr}
  {variant}
  Special point coordinates:
{special_points}
""".format(repr=repr(self),
           variant=self.variant,
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


class Variant:
    variant_desc = """\
Variant name: {name}
  Special point names: {special_point_names}
  Default path: {special_path}
"""

    def __init__(self, name, special_point_names, special_path,
                 special_points=None):
        self.name = name
        self.special_point_names = special_point_names
        self.special_path = special_path
        if special_points is not None:
            special_points = dict(special_points)
            for key, value in special_points.items():
                special_points[key] = np.array(value)
        self.special_points = special_points
        # XXX Should make special_points available as a single array as well
        # (easier to transform that way)

    def __str__(self):
        return self.variant_desc.format(**vars(self))


bravais_names = []
bravais_lattices = {}


def bravaisclass(longname, parameters, variants):
    """Decorator for Bravais lattice classes.

    This sets a number of class variables and processes the information
    about different variants into a list of Variant objects."""

    def decorate(cls):
        btype = cls.__name__
        cls.name = btype
        cls.longname = longname
        cls.parameters = tuple(parameters)
        cls.variant_names = []
        cls.variants = {}

        for [name, special_point_names, special_path,
             special_points] in variants:
            cls.variant_names.append(name)
            cls.variants[name] = Variant(name, special_point_names,
                                         special_path, special_points)

        # Register in global list and dictionary
        bravais_names.append(btype)
        bravais_lattices[btype] = cls
        return cls

    return decorate


class UnconventionalLattice(ValueError):
    pass


class Cubic(BravaisLattice):
    """Abstract class for cubic lattices."""
    def __init__(self, a):
        BravaisLattice.__init__(self, a=a)

@bravaisclass('cubic', 'a',
              [['CUB', 'GXRM', 'GXMGRX,MR', ibz_points['cubic']]])
class CUB(Cubic):
    def _cell(self, a):
        return a * np.eye(3)

@bravaisclass('face-centered cubic', 'a',
              [['FCC', 'GKLUWX', 'GXWKGLUWLK,UX', ibz_points['fcc']]])
class FCC(Cubic):
    def _cell(self, a):
        return 0.5 * np.array([[0., a, a], [a, 0, a], [a, a, 0]])

@bravaisclass('body-centered cubic', 'a',
              [['BCC', 'GHPN', 'GHNGPH,PN', ibz_points['bcc']]])
class BCC(Cubic):
    def _cell(self, a):
        return 0.5 * np.array([[-a, a, a], [a, -a, a], [a, a, -a]])

@bravaisclass('tetragonal', 'ac',
              [['TET', 'GAMRXZ', 'GXMGZRAZ,XR,MA',
                ibz_points['tetragonal']]])
class TET(BravaisLattice):
    def __init__(self, a, c):
        BravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        return np.diag(np.array([a, a, c]))

# XXX in BCT2 we use S for Sigma.
# Also in other places I think
@bravaisclass('body-centered tetragonal', 'ac',
              [['BCT1', 'GMNPXZZ1', 'GXMGZPNZ1M,XP', None],
               ['BCT2', 'GNPSS1XYY1Z', 'GXYSGZS1NPY1Z,XP', None]])
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
        check_orc(a, b, c)
        BravaisLattice.__init__(self, a=a, b=b, c=c)


@bravaisclass('orthorhombic', 'abc',
              [['ORC', 'GRSTUXYZ', 'GXSYGZURTZ,YT,UX,SR',
                ibz_points['orthorhombic']]])
class ORC(Orthorhombic):
    def _cell(self, a, b, c):
        return np.diag([a, b, c]).astype(float)


@bravaisclass('face-centered orthorhombic', 'abc',
              [['ORCF1', 'GAA1LTXX1YZ', 'GYTZGXA1Y,TX1,XAZ,LG', None],
               ['ORCF2', 'GCC1DD1LHH1XYZ', 'GYCDXGZD1HC,C1Z,XH1,HY,LG', None],
               ['ORCF3', 'GAA1LTXX1YZ', 'GYTZGXA1Y,XAZ,LG', None]])
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
        diff = 1.0 / (a * a) - 1.0 / (b * b) - 1.0 / (c * c)
        if abs(diff) < self._eps:
            return 'ORCF3'
        return 'ORCF1' if diff > 0 else 'ORCF2'

@bravaisclass('body-centered orthorhombic', 'abc',
              [['ORCI', 'GLL1L2RSTWXX1YY1Z', 'GXLTWRX1ZGYSW,L1Y,Y1Z', None]])
class ORCI(Orthorhombic):
    def _cell(self, a, b, c):
        return 0.5 * np.array([[-a, b, c], [a, -b, c], [a, b, -c]])

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


@bravaisclass('c-centered orthorhombic', 'abc',
              [['ORCC', 'GAA1RSTXX1YZ', 'GXSRAZGYX1A1TY,ZT', None]])
class ORCC(BravaisLattice):
    def __init__(self, a, b, c):
        # ORCC is the only ORCx lattice with a<b and not a<b<c
        if a >= b:
            raise UnconventionalLattice('Expected a < b, got {}, {}'
                                        .format(a, b, c))
        BravaisLattice.__init__(self, a=a, b=b, c=c)

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


ibz_hex = ibz_points['hexagonal'].copy()
ibz_hex['K'][0] *= -1  # XXX
ibz_hex['H'][0] *= -1

@bravaisclass('hexagonal', 'ac',
              [['HEX', 'GMKALH', 'GMKGALHA,LM,KH', ibz_hex]])
class HEX(BravaisLattice):
    def __init__(self, a, c):
        BravaisLattice.__init__(self, a=a, c=c)

    def _cell(self, a, c):
        x = 0.5 * np.sqrt(3)
        return np.array([[0.5 * a, -x * a, 0], [0.5 * a, x * a, 0],
                         [0., 0., c]])


@bravaisclass('rhombohedral', ('a', 'alpha'),
              [['RHL1', 'GBB1FLL1PP1P2QXZ', 'GLB1,BZGX,QFP1Z,LP', None],
               ['RHL2', 'GFLPP1QQ1Z', 'GPZQGFP1Q1LZ', None]])
class RHL(BravaisLattice):
    def __init__(self, a, alpha):
        if alpha >= 120:
            raise ValueError('Need alpha < 120 degrees, got {}'
                             .format(alpha))
        BravaisLattice.__init__(self, a=a, alpha=alpha)

    def _cell(self, a, alpha):
        alpha *= np.pi / 180
        acosa = a * np.cos(alpha)
        acosa2 = a * np.cos(0.5 * alpha)
        asina2 = a * np.sin(0.5 * alpha)
        acosfrac = acosa / acosa2
        xx = (1 - acosfrac**2)
        assert xx > 0.0
        return np.array([[acosa2, -asina2, 0], [acosa2, asina2, 0],
                         [a * acosfrac, 0, a * xx**0.5]])

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
        raise UnconventionalLattice('Expected a <= c, b <= c, alpha < 90; '
                                    'got a={}, b={}, c={}, alpha={}'
                                    .format(a, b, c, alpha))


@bravaisclass('monoclinic', ('a', 'b', 'c', 'alpha'),
              [['MCL', 'GACDD1EHH1H2MM1M2XYY1Z', 'GYHCEM1AXH1,MDZ,YD', None]])
class MCL(BravaisLattice):
    def __init__(self, a, b, c, alpha):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha)

    def _cell(self, a, b, c, alpha):
        alpha *= np.pi / 180
        return np.array([[a, 0, 0], [0, b, 0],
                         [0, c * np.cos(alpha), c * np.sin(alpha)]])

    def _special_points(self, a, b, c, alpha, variant):
        cosa = np.cos(alpha * _degrees)
        eta = (1 - b * cosa / c) / (2 * np.sin(alpha * _degrees)**2)
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
        return 'MCL'


@bravaisclass('c-centered monoclinic', ('a', 'b', 'c', 'alpha'),
              [['MCLC1', 'GNN1FF1F2F3II1LMXX1X2YY1Z',
                'GYFLI,I1ZF1,YX1,XGN,MG', None],
               ['MCLC2', 'GNN1FF1F2F3II1LMXX1X2YY1Z',
                'GYFLI,I1ZF1,NGM', None],
               ['MCLC3', 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z',
                'GYFHZIF1,H1Y1XGN,MG', None],
               ['MCLC4', 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z',
                'GYFHZI,H1Y1XGN,MG', None],
               ['MCLC5', 'GFF1F2HH1H2II1LMNN1XYY1Y2Y3Z',
                'GYFLI,I1ZHF1,H1Y1XGN,MG', None]])
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


# XXX labels, paths, are all the same.
@bravaisclass('triclinic', ('a', 'b', 'c', 'alpha', 'beta', 'gamma'),
              [['TRI1a', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG', None],
               ['TRI2a', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG', None],
               ['TRI1b', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG', None],
               ['TRI2b', 'GLMNRXYZ', 'XGY,LGZ,NGM,RG', None]])
class TRI(BravaisLattice):
    def __init__(self, a, b, c, alpha, beta, gamma):
        BravaisLattice.__init__(self, a=a, b=b, c=c, alpha=alpha, beta=beta,
                                gamma=gamma)

    def _cell(self, a, b, c, alpha, beta, gamma):
        alpha, beta, gamma = np.array([alpha, beta, gamma])
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


from ase.utils import experimental
@experimental
def get_bravais_lattice(uc, eps=2e-4, _niggli_reduce=False):
    # XXX we need to figure out a way to reduce the cell to standard
    # forms, Niggli or otherwise.  Else we will not always get the
    # right type.

    # orig_uc = uc
    if _niggli_reduce:
        uc, niggli_op = uc.niggli_reduce()

    #uc2 = uc.niggli_reduce()
    #if 0: #np.abs(uc2 - uc).max() > eps:
    #    raise ValueError('Can only get recognize Bravais lattice of '
    #                     'Niggli-reduced cell.')
    #if np.linalg.det(uc.array) < 0:
    #    raise ValueError('Cell should be right-handed')

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
        #for iperm in range(3):
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

    def noduplicates(a):
        sorted_a = np.sort(a)
        smallest_deviation = abs(sorted_a[1:] - sorted_a[:-1]).min()
        return smallest_deviation > eps

    if all_lengths_equal:
        if allclose(angles, 90):
            return CUB(A), {}  #check(CUB, A)
        if allclose(angles, 60):
            return FCC(np.sqrt(2) * A), {}
            #return check(FCC, np.sqrt(2) * A)
        if allclose(angles, np.arccos(-1 / 3) * 180 / np.pi):
            return BCC(2.0 * A / np.sqrt(3)), {}
            #return check(BCC, 2.0 * A / np.sqrt(3))

    if all_lengths_equal and unequal_angle_dir is not None:
        # x = BC_CA_AB[unequal_angle_dir]
        y = BC_CA_AB[(unequal_angle_dir + 1) % 3]

        c = 2.0 * np.sqrt(-y)
        a = np.sqrt(2.0 * A**2 - 0.5 * c**2)
        obj = check(BCT, a, c, axis=-unequal_angle_dir + 2)
        if obj:
            bct = BCT(a, c)
            return bct, {}  # permutation
            #return obj

    if (unequal_angle_dir is not None
          and abs(angles[unequal_angle_dir] - 120) < eps
          and abs(angles[unequal_angle_dir - 1] - 90) < eps):
        a2 = -2 * BC_CA_AB[unequal_scalarprod_dir]
        c = ABC[unequal_scalarprod_dir]
        assert a2 > 0
        lat = HEX(np.sqrt(a2), c)
        return lat, {}  # permutation
        #xxxxxxxx
        #return check(HEX, np.sqrt(a2), c, axis=-unequal_scalarprod_dir + 2)

    if allclose(angles, 90) and unequal_length_dir is not None:
        a = ABC[unequal_length_dir - 1]
        c = ABC[unequal_length_dir]
        xx = check(TET, a, c, axis=-unequal_length_dir + 2)
        assert xx
        return TET(a, c), {}  # permutation

    if (unequal_length_dir is not None
        and abs(BC_CA_AB[unequal_length_dir]) > eps
        and sum(abs(BC_CA_AB) > eps) == 1):

        cdir = unequal_length_dir

        orcc_c = ABC[cdir]
        Y = BC_CA_AB[cdir]
        X = ABC[(cdir + 1) % 3]**2

        #orcc_c = np.sqrt(BC_CA_AB[cdir])
        #sum(abs(BC_CA_AB) > eps)
        #assert Y < 0  # ?
        orcc_a = np.sqrt(2 * (X + Y))
        orcc_b = np.sqrt(2 * (X - Y))
        if orcc_a > orcc_b:
            orcc_a, orcc_b = orcc_b, orcc_a
            # permutation?
        lat = ORCC(orcc_a, orcc_b, orcc_c)
        return lat, {}  # permutation
        #lat.tocell()
        #obj = check(ORCC, a, b, c, axis=2 - unequal_length_dir)
        #xxxxxxxxxxxxxx
        #if obj:
        #    return obj

    if allclose(angles, 90) and all_lengths_different:
        permutation = np.argsort(ABC)
        lat = ORC(*ABC[permutation])
        return lat, {}  # permutation

    if all_lengths_different and noduplicates(BC_CA_AB):
        permutation = np.argsort(BC_CA_AB)
        if all(BC_CA_AB > 0):
            lat = ORCF(*(2 * np.sqrt(BC_CA_AB[permutation])))
            par1 = lat.cellpar()
            par2 = uc.new(uc[permutation]).cellpar()
            if allclose(par1, par2):
                return lat, {}  # XXX permutation

    if all_lengths_equal:
        dims2 = -2 * np.array([BC_CA_AB[1] + BC_CA_AB[2],
                               BC_CA_AB[2] + BC_CA_AB[0],
                               BC_CA_AB[0] + BC_CA_AB[1]])
        if all(dims2 > 0) and noduplicates(dims2):#mindifference > eps:
            dims = np.sqrt(dims2)
            permutation = np.argsort(dims)
            lat = ORCI(*dims[permutation])
            assert allclose(lat.tocell(), uc[permutation])
            return lat, {}  # perm

    if all_lengths_equal:
        cosa = BC_CA_AB[0] / A**2
        rhl_alpha = np.arccos(cosa) * 180 / np.pi
        lat = RHL(A, rhl_alpha)
        return lat, {}
        #obj = check(RHL, A, alpha)
        #xxxxxxxxxxx
        #if obj:
        #    return obj

    #if all_lengths_different and unequal_scalarprod_dir is not None:
    if unequal_scalarprod_dir is not None:
        mcl_alpha = angles[unequal_scalarprod_dir]
        other_two_angles = [angles[unequal_scalarprod_dir - 1],
                            angles[unequal_scalarprod_dir - 2]]
        if allclose(other_two_angles, 90):
            assert mcl_alpha < 90, mcl_alpha
            #cdir = np.argmax(abc)

            #mcl_abc = ABC[np.arange(-3, 0) + unequal_scalarprod_dir]
            a_dir = unequal_scalarprod_dir
            b_dir = (a_dir - 1) % 3
            c_dir = (b_dir - 1) % 3
            mcl_a = ABC[a_dir]
            mcl_b = ABC[b_dir]
            mcl_c = ABC[c_dir]
            if mcl_c < mcl_b:
                b_dir, c_dir = c_dir, b_dir
                mcl_c, mcl_b = mcl_b, mcl_c

            assert mcl_a < mcl_c
            permutation = np.array([a_dir, b_dir, c_dir])
            # Or reverse somehow?

            lat = MCL(mcl_a, mcl_b, mcl_c, mcl_alpha)
            return lat, {}  # permutation
            #cdir = 1 + np.argmax(mcl_abc[1:])
            #mcl_c = mcl_abc[cdir]
            #mcl_b = 

            #mcl_a = 
            #mcl_bc = ABC[unequal_scalarprod_dir - 1]
            #a1 = abc[unequal_scalarprod_dir - 1]
            #a1 = abc[unequal_scalarprod_dir - 2]

            #abc = ABC[np.arange(-3, 0) + unequal_scalarprod_dir]
            #MCL(*abc, alpha=alpha)
        #obj = check(MCL, *abc, alpha=alpha, axis=-unequal_scalarprod_dir)
        #xxxxxxxxxxxxxx
        #if obj:
        #    return obj

    if unequal_length_dir is not None:
        cdir = unequal_length_dir
        mclc_c = ABC[cdir]
        L = ABC[cdir - 1]
        mclc_b = np.sqrt(2 * (L**2 + BC_CA_AB[unequal_length_dir]))
        mclc_a = np.sqrt(4 * L**2 - mclc_b**2)
        cosa = 2 * BC_CA_AB[unequal_length_dir - 1] / (mclc_b * mclc_c)
        mclc_alpha = np.arccos(cosa) * 180 / np.pi
        assert mclc_alpha < 90
        lat = MCLC(mclc_a, mclc_b, mclc_c, mclc_alpha)
        # XXX many things wrong currently
        return lat, {}
        #obj = check(MCLC, a, b, c, alpha, axis=-unequal_length_dir + 2)
        #xxxxxxxxxxx
        #if obj:
        #    return obj

    obj = check(TRI, A, B, C, *angles)
    if obj:
        cls, op = obj
        #xxxxxxxxxxx
        # Should always be true
        return cls(A, B, C, *angles), op

    raise RuntimeError('Cannot recognize cell at all somehow!')


def all_variants():
    """For testing and examples; yield all variants of all lattices."""
    a, b, c = 3., 4., 5.
    alpha = 55.0
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

    rhl1 = RHL(a, alpha=55.0)
    assert rhl1.variant.name == 'RHL1'
    yield rhl1

    rhl2 = RHL(a, alpha=105.0)
    assert rhl2.variant.name == 'RHL2'
    yield rhl2

    yield MCL(a, b, c, alpha)

    mclc1 = MCLC(a, b, c, 80)
    assert mclc1.variant.name == 'MCLC1'
    yield mclc1
    # mclc2 has same special points as mclc1

    mclc3 = MCLC(1.8 * a, b, c * 2, 80)
    assert mclc3.variant.name == 'MCLC3'
    yield mclc3
    # mclc4 has same special points as mclc3

    mclc5 = MCLC(b, b, 1.1 * b, 50)
    assert mclc5.variant.name == 'MCLC5'
    yield mclc5

    def get_tri(kcellpar):
        # We build the TRI lattices from cellpars of reciprocal cell
        icell = Cell.fromcellpar(kcellpar)
        cellpar = Cell(4 * icell.reciprocal()).cellpar()
        return TRI(*cellpar)

    tri1a = get_tri([1., 1.2, 1.4, 120., 110., 100.])
    assert tri1a.variant.name == 'TRI1a'
    yield tri1a

    tri1b = get_tri([1., 1.2, 1.4, 50., 60., 70.])
    assert tri1b.variant.name == 'TRI1b'
    yield tri1b

    tri2a = get_tri([1., 1.2, 1.4, 120., 110., 90.])
    assert tri2a.variant.name == 'TRI2a'
    yield tri2a
    tri2b = get_tri([1., 1.2, 1.4, 50., 60., 90.])
    assert tri2b.variant.name == 'TRI2b'
    yield tri2b
