def test_cfg():
    import numpy as np

    from ase_koopmans.build import molecule
    from ase_koopmans.io import read, write

    a = molecule('CO2')
    f = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    a.set_array('test', f)

    write('test.cfg', a)

    b = read('test.cfg')
    assert np.all(b.get_array('test') == f)

    a.set_momenta(2 * f)
    write('test.cfg', a)

    b = read('test.cfg')
    assert np.all(np.abs(a.get_momenta() - b.get_momenta()) < 1e-3)
