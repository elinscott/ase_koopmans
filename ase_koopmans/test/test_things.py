def test_things():
    from ase_koopmans.dft.kpoints import monkhorst_pack
    from ase_koopmans.units import Hartree, Bohr, kJ, mol, kcal, kB, fs
    from ase_koopmans.build import bulk

    assert [0, 0, 0] in  monkhorst_pack((1, 3, 5)).tolist()
    assert [0, 0, 0] not in  monkhorst_pack((1, 3, 6)).tolist()
    assert len(monkhorst_pack((3, 4, 6))) == 3 * 4 * 6

    print(Hartree, Bohr, kJ/mol, kcal/mol, kB * 300, fs, 1 / fs)

    hcp = bulk('X', 'hcp', a=1) * (2, 2, 1)
    assert abs(hcp.get_distance(0, 3, mic=True) - 1) < 1e-12
    assert abs(hcp.get_distance(0, 4, mic=True) - 1) < 1e-12
    assert abs(hcp.get_distance(2, 5, mic=True) - 1) < 1e-12
