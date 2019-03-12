import numpy as np
from ase.geometry.cell import Cell
from ase.geometry.bravais import bravais_lattices as bravais1
from ase.build import bulk

bravais = {}
for name in bravais1:
    bravais[name.lower()] = bravais1[name]

def check_single(name, cell):
    c = Cell(cell)
    lattice, _ = c.bravais()
    name1 = lattice.name.lower()
    ok = name.split('@')[0] == name1
    print(name, '-->', name1, 'OK' if ok else 'ERR', c.cellpar())
    assert ok

def check(name, cell):
    cell = Cell(cell).array
    #cell = Cell(cell).array
    # Check all three positive permutations:
    check_single(name + '@012', cell[[0, 1, 2]])
    # check_single(name + '@201', cell[[2, 0, 1]])
    # check_single(name + '@120', cell[[1, 2, 0]])

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
# check('orc', bravais['orc'](4., 5., 3.).tocell())
check('orcf', bravais['orcf'](4., 5., 7.).tocell())
check('orci', bravais['orci'](2., 5., 6.).tocell())
check('orcc', bravais['orcc'](3., 4., 5.).tocell())
check('hex', bravais['hex'](5., 6.).tocell())
check('rhl', bravais['rhl'](4., 54.).tocell())
check('mcl', bravais['mcl'](2., 3., 4., 62.).tocell())
check('mclc', bravais['mclc'](3., 4., 5., 70.).tocell())
check('tri', bravais['tri'](7., 6., 5., 65., 70., 80.).tocell())
