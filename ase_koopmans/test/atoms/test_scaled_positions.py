def test_scaled_positions():
    from ase_koopmans import Atoms
    assert Atoms('X', [(-1e-35, 0, 0)], pbc=True).get_scaled_positions()[0, 0] < 1
