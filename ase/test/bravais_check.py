import numpy as np
from ase.geometry.cell import Cell
from ase.geometry.bravais import bravais_lattices
from ase.build import bulk, fcc111

bravais = {}
for name in bravais_lattices:
    bravais[name.lower()] = bravais_lattices[name]

def check_single(name, cell, pbc=None):
    c = Cell(cell, pbc=pbc)
    try:
        lattice = c.get_bravais_lattice()
    except RuntimeError:
        print('error checking {}'.format(name))
        raise
    name1 = lattice.name.lower()
    latname = name.split('@')[0]
    ok = latname == name1
    print(name, '-->', name1, 'OK' if ok else 'ERR', c.cellpar())
    assert ok, 'Expected {latname} but found {name1}'.format(latname, name1)


def check(name, cell):
    if isinstance(cell, Cell):
        pbc = cell.pbc
    else:
        pbc = np.array([True] * 3)
    cell = Cell(cell).array

    # Check all three positive permutations:
    check_single(name + '@012', cell[[0, 1, 2]], pbc=pbc)
    # XXX 2D lattice determination not implemented for arbitrary permutations
    if pbc.sum() == 3:
        check_single(name + '@201', cell[[2, 0, 1]], pbc=pbc)
        check_single(name + '@120', cell[[1, 2, 0]], pbc=pbc)


check('cub', bravais['cub'](3.3).tocell())
check('fcc', bravais['fcc'](3.4).tocell())
check('fcc', bulk('Au').cell)
check('bcc', bravais['bcc'](3.5).tocell())
check('bcc', bulk('Fe').cell)
check('tet', bravais['tet'](4., 5.).tocell())
check('tet', np.diag([4., 5., 5.]))
check('tet', np.diag([5., 4., 5.]))
check('tet', np.diag([5., 5., 4.]))
check('bct', bravais['bct'](3., 4.).tocell())
check('orc', bravais['orc'](3., 4., 5.).tocell())
check('orcf', bravais['orcf'](4., 5., 7.).tocell())
check('orci', bravais['orci'](2., 5., 6.).tocell())
check('orcc', bravais['orcc'](3., 4., 5.).tocell())
check('hex', fcc111('Au', size=(1, 1, 3), periodic=True).cell)
check('hex', bravais['hex'](5., 6.).tocell())
check('rhl', bravais['rhl'](4., 54.).tocell())
check('mcl', bravais['mcl'](2., 3., 4., 62.).tocell())
check('mclc', bravais['mclc'](3., 4., 5., 75.).tocell())
#check('tri', bravais['tri'](7., 6., 5., 65., 70., 80.).tocell())

# For 2D materials we have to check both the tocell() method
# but also for realistic cell nonzero nonperiodic axis.
check('sqr', bravais['sqr'](3.).tocell())
check('sqr', Cell(np.diag([3., 3., 10.]),
                  pbc=np.array([True, True, False])))

check('crect', bravais['crect'](3., 40).tocell())

alpha = 40 / 360 * 2 * np.pi
a = 3
x = np.cos(alpha)
y = np.sin(alpha)

crectcell = np.array([[a, 0, 0],
                      [a * x, a * y, 0],
                      [0, 0, 10]])
check('crect', Cell(crectcell, pbc=np.array([True, True, False])))

check('rect', bravais['rect'](3., 4.).tocell())
check('rect', Cell(np.diag([3., 4., 10.]),
                   pbc=np.array([True, True, False])))

check('hex2d', bravais['hex2d'](3.).tocell())
x = 0.5 * np.sqrt(3)
hexcell = np.array([[a, 0, 0],
                    [-0.5 * a, x * a, 0],
                    [0., 0., 0.]])
check('hex2d', Cell(hexcell, pbc=np.array([True, True, False])))

check('obl', bravais['obl'](3., 4., 40).tocell())

b = 4
x = np.cos(alpha)
y = np.sin(alpha)
oblcell = np.array([[a, 0, 0],
                    [b * x, b * y, 0],
                    [0, 0, 10]])
check('obl', Cell(oblcell, pbc=np.array([True, True, False])))
