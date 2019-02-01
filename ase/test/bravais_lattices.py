from ase.geometry.bravais import bravais_lattices, _test_all_variants

for name in bravais_lattices:
    latcls = bravais_lattices[name]
    assert latcls.type == name
    assert latcls.type is not None
    assert latcls.name is not None
    for par in latcls.parameters:
        assert par in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']

    # Check variants?


for lat in _test_all_variants():
    print(lat.variant)
    for par in lat.parameters:
        print(par, getattr(lat, par))

    print('cell', lat.tocell())  # tocell().bravais()[0] should evaluate back to same lattice
    print('cellpar', lat.cellpar())
    print('special path', lat.special_path)
    arr = lat.get_special_points_array()
    assert arr.shape == (len(lat.special_point_names), 3)

    dct = lat.get_special_points()
    assert len(dct) == len(lat.special_point_names)
    print(lat)
