def test_parallel():
    from ase_koopmans import Atoms
    from ase_koopmans.io import read, write
    from ase_koopmans.parallel import world

    n = world.rank + 1
    a = Atoms('H' * n)
    name = 'H{}.xyz'.format(n)
    write(name, a, parallel=False)
    b = read(name, parallel=False)
    assert n == len(b)
