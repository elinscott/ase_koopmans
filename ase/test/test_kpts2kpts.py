from ase.calculators.calculator import kpts2kpts
from ase.lattice import all_variants
from ase import Atoms


# This function tests whether giving a bandpath
# and an atoms object with a completed cell to
# kpts2kpts actually produces a band path with
# the same special points in the end. If this
# isn't fulfilled then it means that something
# has gone wrong (most likely that the bravais
# lattice wasn't correctly identified).
for lat in all_variants():
    print()
    print(lat)
    bandpath = lat.bandpath()
    a = Atoms()
    a.cell = lat.tocell().complete()
    a.pbc = [True if i < lat.ndim else False for i in range(3)]
    path = {'path': bandpath.path}
    bandpath2 = kpts2kpts(path, atoms=a)
    print('cell', a.cell)
    print('Original', bandpath)
    print('path', path)
    print('Produced by kpts2kpts', bandpath2)
    sp = set(bandpath.special_points.keys())
    sp2 = set(bandpath2.special_points.keys())
    intsp = sp & sp2
    msg = ('Input and output bandpath from kpts2kpts dont agree!\n'
           f'Input: {bandpath}\n Output: {bandpath2}')
    assert not sp - intsp, msg
    assert not sp2 - intsp, msg
