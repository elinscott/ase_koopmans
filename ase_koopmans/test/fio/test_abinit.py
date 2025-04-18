def test_abinit():
    import numpy as np

    def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
        """Replacement for old numpy.testing.utils.array_almost_equal."""
        return (np.abs(a1 - a2) < tol).all()

    # this test should be run with abinit!
    from ase_koopmans.calculators.emt import EMT

    from ase_koopmans.io import read, write

    from ase_koopmans.build import molecule

    m1 = molecule('O2', pbc=True)
    m1.center(2.0)

    write('abinit_save.in', images=m1, format='abinit-in')

    m1.calc = EMT()
    e1 = m1.get_potential_energy()
    f1 = m1.get_forces()

    m2 = read('abinit_save.in', format='abinit-in')

    m2.calc = EMT()
    e2 = m2.get_potential_energy()
    f2 = m1.get_forces()

    # assume atoms definitions are the same if energy/forces are the same: can we do better?
    assert abs(e1-e2) < 1.e-6, str(e1) + ' ' + str(e2)
    assert array_almost_equal(f1, f2, tol=1.e-6)
