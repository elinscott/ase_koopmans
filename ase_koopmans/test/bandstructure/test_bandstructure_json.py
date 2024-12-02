def test_bandstructure_json():
    from ase_koopmans.build import bulk
    from ase_koopmans.spectrum.band_structure import calculate_band_structure, BandStructure
    from ase_koopmans.io.jsonio import read_json
    from ase_koopmans.calculators.test import FreeElectrons

    atoms = bulk('Au')
    lat = atoms.cell.get_bravais_lattice()
    path = lat.bandpath(npoints=100)

    atoms.calc = FreeElectrons()

    bs = calculate_band_structure(atoms, path)
    bs.write('bs.json')
    bs.path.write('path.json')

    bs1 = read_json('bs.json')
    bs2 = BandStructure.read('bs.json')
    path1 = read_json('path.json')
    assert type(bs1) == type(bs)  # noqa
    assert type(bs2) == type(bs)  # noqa
    assert type(path1) == type(bs.path)  # noqa
