def test_issue276():
    import warnings

    import numpy as np

    from ase_koopmans.io import read, write
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.build import bulk

    at = bulk("Cu")
    at.rattle()
    at.calc = EMT()
    f = at.get_forces()

    write("tmp.xyz", at)
    at2 = read("tmp.xyz")
    f2 = at.get_forces()

    assert np.abs(f - f2).max() < 1e-6

    with warnings.catch_warnings(record=True) as w:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        write("tmp2.xyz", at2)
        assert len(w) == 2
        assert ('overwriting array' in str(w[0].message))
        assert ('overwriting array' in str(w[1].message))

    at3 = read("tmp2.xyz")
    f3 = at3.get_forces()
    assert np.abs(f - f3).max() < 1e-6
