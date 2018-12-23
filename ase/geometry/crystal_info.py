import numpy as np
#from ase.geometry.cell import Cell#, get_bravais_lattice
#from ase.build.tools import niggli_reduce_cell

special_paths = {
    'cubic': 'GXMGRX,MR',
    'fcc': 'GXWKGLUWLK,UX',
    'bcc': 'GHNGPH,PN',
    'tetragonal': 'GXMGZRAZXR,MA',
    'orthorhombic': 'GXSYGZURTZ,YT,UX,SR',
    'hexagonal': 'GMKGALHA,LM,KH',
    'monoclinic': 'GYHCEM1AXH1,MDZ,YD',
    'rhombohedral type 1': 'GLB1,BZGX,QFP1Z,LP',
    'rhombohedral type 2': 'GPZQGFP1Q1LZ'}


ibz_points = {'cubic': {'Gamma': [0, 0, 0],
                        'X': [0, 0 / 2, 1 / 2],
                        'R': [1 / 2, 1 / 2, 1 / 2],
                        'M': [0 / 2, 1 / 2, 1 / 2]},
              'fcc': {'Gamma': [0, 0, 0],
                      'X': [1 / 2, 0, 1 / 2],
                      'W': [1 / 2, 1 / 4, 3 / 4],
                      'K': [3 / 8, 3 / 8, 3 / 4],
                      'U': [5 / 8, 1 / 4, 5 / 8],
                      'L': [1 / 2, 1 / 2, 1 / 2]},
              'bcc': {'Gamma': [0, 0, 0],
                      'H': [1 / 2, -1 / 2, 1 / 2],
                      'N': [0, 0, 1 / 2],
                      'P': [1 / 4, 1 / 4, 1 / 4]},
              'hexagonal': {'Gamma': [0, 0, 0],
                            'M': [0, 1 / 2, 0],
                            'K': [-1 / 3, 1 / 3, 0],
                            'A': [0, 0, 1 / 2],
                            'L': [0, 1 / 2, 1 / 2],
                            'H': [-1 / 3, 1 / 3, 1 / 2]},
              'tetragonal': {'Gamma': [0, 0, 0],
                             'X': [1 / 2, 0, 0],
                             'M': [1 / 2, 1 / 2, 0],
                             'Z': [0, 0, 1 / 2],
                             'R': [1 / 2, 0, 1 / 2],
                             'A': [1 / 2, 1 / 2, 1 / 2]},
              'orthorhombic': {'Gamma': [0, 0, 0],
                               'R': [1 / 2, 1 / 2, 1 / 2],
                               'S': [1 / 2, 1 / 2, 0],
                               'T': [0, 1 / 2, 1 / 2],
                               'U': [1 / 2, 0, 1 / 2],
                               'X': [1 / 2, 0, 0],
                               'Y': [0, 1 / 2, 0],
                               'Z': [0, 0, 1 / 2]}}



class CrystalSystem:
    def __init__(self, name, lattices, space_groups):
        self.name = name
        self.lattices = tuple(lattices)
        self.space_groups = space_groups

    def __repr__(self):
        return ('{}(name={},lattices={},space_groups={})'
                .format(self.__class__.__name__, self.lattices,
                        self.space_groups))

def _get_crystal_systems():
    crystal_system_data = [('triclinic', ['tri'], range(1, 3)),
                           ('monoclinic', ['mcl', 'mclc'], range(3, 16)),
                           ('orthorhombic', ['tet', 'bct'], range(16, 75)),
                           ('tetragonal', ['tet', 'bct'], range(75, 142)),
                           ('trigonal', ['rhl', 'hex'], range(143, 168)),
                           ('hexagonal', ['hex'], range(168, 195)),  # D'oh!
                           ('cubic', ['cub', 'fcc', 'bcc'], range(195, 231))]
    from collections import OrderedDict

    d = OrderedDict()
    for name, lattices, spgs in crystal_system_data:
        d[name] = CrystalSystem(name, lattices, spgs)
    return d

crystal_systems = _get_crystal_systems()

# Shortcut for defining functions where we do not have several variants
class _special_points:
    def __init__(self, name, labels):
        self.name = name
        self.labels = labels

    def __call__(self):
        return (None,  # variant; in these cases we do not have a variant
                self.labels, special_paths[self.name].replace(',', ' '),
                [ibz_points[self.name][label] for label in self.labels])


cub_special_points = _special_points('cubic', 'GXRM')
fcc_special_points = _special_points('fcc', 'GKLUWX')
bcc_special_points = _special_points('bcc', 'GHPN')
tet_special_points = _special_points('tetragonal', 'GAMRXZ')

def bct_special_points(a, c):
    a2 = a * a
    c2 = c * c
    eta = .25 * (1 + c2 / a2)

    variant = 1 if c < a else 2

    if variant == 1:
        names = 'GMNPXSS1'
        paths = 'GXMGSPNS1M XP'
        points = [[0,0,0],
                  [-.5, .5, .5],
                  [0.,.5,0.],
                  [.25, .25, .25],
                  [0.,0.,.5],
                  [eta,eta,-eta],
                  [-eta,1-eta,eta]]
    else:
        zeta = 0.5 * a2 / c2
        names = 'GNPSS1XYY1Z'
        paths = 'GXYSGZS1NPY1Z XP'
        points = [[0.,.0,0.],
                  [0.,.5,0.],
                  [.25,.25,.25],
                  [-eta,eta,eta],
                  [eta,1-eta,-eta],
                  [0.,0.,.5],
                  [-zeta,zeta,.5],
                  [.5,.5,-zeta],
                  [.5,.5,-.5]]
    return variant, names, paths, points


orc_special_points = _special_points('orthorhombic', 'GRSTUXYZ')

def orcf_special_points(a, b, c):
    a2 = a * a
    b2 = b * b
    c2 = c * c
    xminus = 0.25 * (1 + a2 / b2 - a2 / c2)
    xplus = 0.25 * (1 + a2 / b2 + a2 / c2)

    variant = get_orcf_variant(a, b, c)

    if variant == 1 or variant == 3:
        names = 'GAA1LTXX1YZ'
        paths = 'GYTZGXA1Y TX1 XAZ LG'

        zeta = xminus
        eta = xplus

        # re.findall('[A-Z]\d?', names) <-- convert to list
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

        names = 'GCC1DD1LHH1XYZ'
        paths = 'GYCDXGZD1HC C1Z XH1 HY LG'
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

    return variant, names, paths, points

def orci_special_points(a, b, c):
    a2 = a**2
    b2 = b**2
    c2 = c**2

    zeta = .25 * (1 + a2 / c2)
    eta = .25 * (1 + b2 / c2)
    delta = .25 * (b2 - a2) / c2
    mu = .25 * (a2 + b2) / c2

    names = 'GLL1L2RSTWXX1YY1Z'
    paths = 'GXLTWRX1ZGYSW L1Y Y1Z'
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
    return 1, names, paths, points

def orcc_special_points(a, b, c):
    zeta = .25 * (1 + a * a / (b * b))
    names = 'GAA1RSTXX1YZ'
    paths = 'GXSRAZGYX1A1TY ZT'
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
    return 1, names, paths, points

hex_special_points = _special_points('hexagonal', 'GMKALH')

def rhl_special_points(a, alpha):
    cosa = np.cos(alpha)
    eta = (1 + 4 * cosa) / (2 + 4 * cosa)
    nu = .75 - 0.5 * eta

    variant = 1 if alpha < 90 else 2

    if variant == 1:
        names = 'GBB21FLL1PP1P2QXZ'
        paths = 'GLB1 BZGX QFP1Z LP'
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
        names = 'GFLPP1QQ1Z'
        paths = 'GPZQGFP1Q1LZ'
        points = [[0,0,0],
                  [.5,-.5,0],
                  [.5,0,0],
                  [1-nu,-nu,1-nu],
                  [nu,nu-1,nu-1],
                  [eta,eta,eta],
                  [1-eta,-eta,-eta],
                  [.5,-.5,.5]]
    return variant, names, paths, points

def mcl_special_points(a, b, c, alpha):
    assert a < c
    assert b < c
    assert alpha < 90

    cosa = np.cos(alpha)
    eta = (1 - b * cosa / c) / (2 * np.sin(alpha)**2)
    nu = .5 - eta * c * cosa / b

    names = 'GACDD1EHH1H2MM1M2XYY1Z'
    paths = 'GYHCEM1AXH1 MDZ YD'
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
    return 1, names, paths, points

def mclc_special_points(a, b, c, alpha, eps=2e-4):
    assert a <= c
    assert b <= c
    assert alpha < 90

    a2 = a * a
    b2 = b * b
    # c2 = c * c
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    sina2 = sina**2

    from ase.geometry.cell import mclc
    cell = mclc(a, b, c, alpha)
    lengths_angles = Cell(np.linalg.inv(cell)).cellpar()

    kgamma = lengths_angles[-1]

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

    paths = {1: 'GYFLI I1ZF1 YX1 XGN MG',
             2: 'GYFLI I1ZF1 NGM',
             3: 'GYFHZIF1 H1Y1XGN MG',
             4: 'GYFHZI H1Y1XGN MG',
             5: 'GYFLI I1ZHF1 H1Y1XGN MG'}

    if variant == 1 or variant == 2:
        zeta = (2 - b * cosa / c) / (4 * sina2)
        eta = 0.5 + 2 * zeta * c * cosa / b
        psi = .75 - a2 / (4 * b2 * sina * sina)
        phi = psi + (.75 - psi) * b * cosa / c

        names = 'GNN1FF1F2F3II1LMXX1X2YY1Z'
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

        names = 'GFF1F2HH1H2IMNN1XYY1Y2Y3Z'
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
        names = 'GFF1F2HH1H2II1LMNN1XYY1Y2Y3Z'
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

    return variant, names, paths[variant], points

def tri_special_points(a, b, c, alpha, beta, gamma, eps=2e-4):
    c = Cell.new([a, b, c, alpha, beta, gamma])
    reciprocal = c.reciprocal()
    from ase.geometry.cell import get_cell_lengths_and_angles
    ka, kb, kc, kalpha, kbeta, kgamma = get_cell_lengths_and_angles(reciprocal)
    # lengths = np.array([ka, kb, kc])
    angles = np.array([kalpha, kbeta, kgamma])

    paths = {'1a': 'XGY LGZ NGM RG',
             '2a': 'XGY LGZ NGM RG',
             '1b': 'XGY LGZ NGM RG',
             '2b': 'XGY LGZ NGM RG'}  # Eh, they are the same!

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

    names = 'GLMNRXYZ'
    if var == '1a' or var == '2a':
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

    return var, names, paths[var], points

def _get_special_points():
    # name/variant | special point names | special paths | special points
    d = [('cub', 'GMRX', ['GXMGRX', 'MR'],
          [[0, 0, 0],
           [1 / 2, 1 / 2, 0],
           [1 / 2, 1 / 2, 1 / 2],
           [0, 1 / 2, 0]]),
         ('fcc', 'GKLUWX', ['GHNGPH', 'PN'],
          [[0, 0, 0],
           [3 / 8, 3 / 8, 3 / 4],
           [1 / 2, 1 / 2, 1 / 2],
           [5 / 8, 1 / 4, 5 / 8],
           [1 / 2, 1 / 4, 3 / 4],
           [1 / 2, 0, 1 / 2]]),
          ('bcc'),
          ('tet'),
          ('orc'),
          ('orcf1', 'GAA1LTXX1YZ', ['GYTZGXA1Y', 'TX1', 'XAZ', 'LG'],
          ),
         # orcf3: same as orcf1
          ('orcc'),
          (),
          (),
          (),
          (),
          (),
          ('mcl'),
          ('mclc'),
          ('tri'),
          ()]
    return d


def _get_bravais_data():
    def get_orcf_variant(a, b, c, eps=2e-4):
        diff = 1.0 / (a * a) - 1.0 / (b * b) - 1.0 / (c * c)
        if abs(diff) < eps:
            return '3'
        return '1' if diff > 0 else '2'

    def get_mcl_variant(a, b, c, alpha):
        raise NotImplementedError  # grumble depends on kgamma

    def get_tri_variant(a, b, c, alpha, beta, gamma):
        raise NotImplementedError

    # Short name | long name | number of variant | get_variant (callable)
    d = [('cub', 'cubic', 1, None),
         ('fcc', 'face-centered cubic', 1, None),
         ('bcc', 'body-centered cubic', 1, None),
         ('tet', 'tetragonal', 2, lambda a, c, eps=None: '2' if c > a else '1'),
         ('bct', 'body-centered tetragonal', 1, None),
         ('orc', 'orthorhombic', 1, None),
         ('orcf', 'face-centered orthorhombic', 3, get_orcf_variant),
         ('orci', 'body-centered orthorhombic', 1, None),
         ('orcc', 'c-centered orthorhombic', 1, None),
         ('hex', 'hexagonal', 1, None),
         ('rhl', 'rhombohedral', 2,
          lambda a, alpha: '1' if alpha < 90.0 else '2'),
         ('mcl', 'monoclinic', 1, None),
         ('mclc', 'c-centered monoclinic', 5, get_mcl_variant),
         ('tri', 'triclinic', 4, get_tri_variant)]
    return d


class ExtendedCellInfo:
    def __init__(self, **kwargs):
        self._d = kwargs

    @property
    def bravais(self):  # could be called 'type'
        return self._d['bravais']

    @property
    def parameters(self):
        return self._d['parameters'].copy()

    def __repr__(self):
        kwstring = ', '.join('{}={}'.format(k, v) for k, v
                             in sorted(self._d.items()))
        return '{}({})'.format(self.__class__.__name__, kwstring)

    @property
    def cell(self):
        return self.bravais(**self.parameters)

    #def bandpath(path=None, npoints):
    #    if path is None:
    #        path = self.special_bandpaths[0]

    #    return BandPath(...)

    @property
    def special_bandpaths(self):
        return ...

    @property
    def special_points(self):
        return ...

    @property
    def variant(self):
        return ...

    # bandpath()
    # default bandpaths
    # special points
    # bravais variant


def analyse_cell(cell, eps=2e-4):
    rcell, _ = niggli_reduce_cell(cell)
    bravais, latticepar = get_bravais_lattice(Cell(rcell), eps=eps)

    print(bravais)
    print(latticepar)
    info = ExtendedCellInfo(bravais=bravais, parameters=latticepar,
                            pbc=cell.pbc)
    return info
