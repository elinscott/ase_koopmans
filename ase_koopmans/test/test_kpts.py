def test_kpts():
    from ase_koopmans.dft.kpoints import bandpath
    import numpy as np
    print(bandpath('GX,GX', np.eye(3), 6))
