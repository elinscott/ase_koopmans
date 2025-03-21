def test_castep_interface():
    """Simple shallow test of the CASTEP interface"""
    import os
    import re
    import tempfile
    import warnings

    import numpy as np
    import ase_koopmans
    import ase_koopmans.lattice.cubic
    from ase_koopmans.calculators.castep import (Castep, CastepOption,
                                        CastepParam, CastepCell,
                                        make_cell_dict, make_param_dict,
                                        CastepKeywords,
                                        create_castep_keywords,
                                        import_castep_keywords,
                                        CastepVersionError)

    # XXX on porting this test to pytest it wasn't skipped as it should be.
    # At any rate it failed then.  Maybe someone should look into that ...
    #
    # Hence, call the constructor to trigger our test skipping hack:
    Castep()

    tmp_dir = tempfile.mkdtemp()

    # We have fundamentally two sets of tests: one if CASTEP is present, the other
    # if it isn't
    has_castep = False
    # Try creating and importing the castep keywords first
    try:
        create_castep_keywords(
            castep_command=os.environ['CASTEP_COMMAND'],
            path=tmp_dir,
            fetch_only=20)
        has_castep = True  # If it worked, it must be present
    except KeyError:
        print('Could not find the CASTEP_COMMAND environment variable - please_koopmans'
              ' set it to run the full set of Castep tests')
    except CastepVersionError:
        print('Invalid CASTEP_COMMAND provided - please_koopmans set the correct one to '
              'run the full set of Castep tests')

    try:
        castep_keywords = import_castep_keywords(
            castep_command=os.environ.get('CASTEP_COMMAND', ''))
    except CastepVersionError:
        castep_keywords = None

    # Start by testing the fundamental parts of a CastepCell/CastepParam object
    boolOpt = CastepOption('test_bool', 'basic', 'defined')
    boolOpt.value = 'TRUE'
    assert boolOpt.raw_value is True

    float3Opt = CastepOption('test_float3', 'basic', 'real vector')
    float3Opt.value = '1.0 2.0 3.0'
    assert np.isclose(float3Opt.raw_value, [1, 2, 3]).all()

    # Generate a mock keywords object
    mock_castep_keywords = CastepKeywords(make_param_dict(), make_cell_dict(),
                                          [], [], 0)
    mock_cparam = CastepParam(mock_castep_keywords, keyword_tolerance=2)
    mock_ccell = CastepCell(mock_castep_keywords, keyword_tolerance=2)

    # Test special parsers
    mock_cparam.continuation = 'default'
    mock_cparam.reuse = 'default'
    assert mock_cparam.reuse.value is None

    mock_ccell.species_pot = ('Si', 'Si.usp')
    mock_ccell.species_pot = ('C',  'C.usp')
    assert 'Si Si.usp' in mock_ccell.species_pot.value
    assert 'C C.usp' in mock_ccell.species_pot.value
    symops = (np.eye(3)[None], np.zeros(3)[None])
    mock_ccell.symmetry_ops = symops
    assert """1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
    0.0 0.0 0.0""" in mock_ccell.symmetry_ops.value

    # check if the CastepOpt, CastepCell comparison mechanism works
    if castep_keywords:
        p1 = CastepParam(castep_keywords)
        p2 = CastepParam(castep_keywords)

        assert p1._options == p2._options

        p1._options['xc_functional'].value = 'PBE'
        p1.xc_functional = 'PBE'

        assert p1._options != p2._options

    c = Castep(directory=tmp_dir, label='test_label', keyword_tolerance=2)
    if castep_keywords:
        c.xc_functional = 'PBE'
    else:
        c.param.xc_functional = 'PBE'  # In "forgiving" mode, we need to specify

    lattice = ase_koopmans.lattice.cubic.BodyCenteredCubic('Li')

    print('For the sake of evaluating this test, warnings')
    print('about auto-generating pseudo-potentials are')
    print('normal behavior and can be safely ignored')

    lattice.calc = c

    param_fn = os.path.join(tmp_dir, 'myParam.param')
    param = open(param_fn, 'w')
    param.write('XC_FUNCTIONAL : PBE #comment\n')
    param.write('XC_FUNCTIONAL : PBE #comment\n')
    param.write('#comment\n')
    param.write('CUT_OFF_ENERGY : 450.\n')
    param.close()
    c.merge_param(param_fn)

    assert c.calculation_required(lattice)
    if has_castep:
        assert c.dryrun_ok()

    c.prepare_input_files(lattice)

    # detecting pseudopotentials tests

    # typical filenames
    files = ['Ag_00PBE.usp',
             'Ag_00.recpot',
             'Ag_C18_PBE_OTF.usp',
             'ag-optgga1.recpot',
             'Ag_OTF.usp',
             'ag_pbe_v1.4.uspp.F.UPF',
             'Ni_OTF.usp',
             'fe_pbe_v1.5.uspp.F.UPF',
             'Cu_01.recpot']

    pp_path = os.path.join(tmp_dir, 'test_pp')
    os.makedirs(pp_path)

    for f in files:
        with open(os.path.join(pp_path, f), 'w') as _f:
            _f.write('DUMMY PP')


    c = Castep(directory=tmp_dir, label='test_label_pspots',
               castep_pp_path=pp_path)
    c._pedantic=True
    atoms=ase_koopmans.build.bulk('Ag')
    atoms.calc = c

    # I know, unittest would be nicer... maybe at a later point

    # disabled, but may be useful still
    # try:
    # # this should yield no files
    # atoms.calc.find_pspots(suffix='uspp')
    # raise AssertionError
    #    # this should yield no files
    #    atoms.calc.find_pspots(suffix='uspp')
    #    raise AssertionError
    # except RuntimeError as e:
    # #print(e)
    # pass
    #     # print(e)
    #     pass

    try:
        # this should yield non-unique files
        atoms.calc.find_pspots(suffix = 'recpot')
        raise AssertionError
    except RuntimeError:
        pass

    # now let's see if we find all...
    atoms.calc.find_pspots(pspot = '00PBE', suffix = 'usp')
    assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_00PBE.usp'

    atoms.calc.find_pspots(pspot = '00', suffix = 'recpot')
    assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_00.recpot'

    atoms.calc.find_pspots(pspot = 'C18_PBE_OTF', suffix = 'usp')
    assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_C18_PBE_OTF.usp'

    atoms.calc.find_pspots(pspot = 'optgga1', suffix = 'recpot')
    assert atoms.calc.cell.species_pot.value.split()[-1] == 'ag-optgga1.recpot'

    atoms.calc.find_pspots(pspot = 'OTF', suffix = 'usp')
    assert atoms.calc.cell.species_pot.value.split()[-1] == 'Ag_OTF.usp'

    atoms.calc.find_pspots(suffix = 'UPF')
    assert (atoms.calc.cell.species_pot.value.split()[-1] ==
            'ag_pbe_v1.4.uspp.F.UPF')


    # testing regular workflow
    c = Castep(directory=tmp_dir, label='test_label_pspots',
               castep_pp_path=pp_path, find_pspots=True, keyword_tolerance=2)
    c._build_missing_pspots = False
    atoms = ase_koopmans.build.bulk('Ag')
    atoms.calc = c

    # this should raise an error due to ambuiguity
    try:
        c._fetch_pspots()
        raise AssertionError
    except RuntimeError:
        pass

    for e in ['Ni', 'Fe', 'Cu']:
        atoms = ase_koopmans.build.bulk(e)
        atoms.calc = c
        c._fetch_pspots()

    # test writing to file
    tmp_dir = os.path.join(tmp_dir, 'input_files')
    c = Castep(directory=tmp_dir,
               find_pspots=True, castep_pp_path=pp_path, keyword_tolerance=2)
    c._label = 'test'
    atoms = ase_koopmans.build.bulk('Cu')
    atoms.calc = c
    c.prepare_input_files()

    with open(os.path.join(tmp_dir, 'test.cell'), 'r') as f:
        assert re.search(r'Cu Cu_01\.recpot', ''.join(f.readlines())) is not None


    # test keyword conflict management
    c = Castep(cut_off_energy=300.)
    assert float(c.param.cut_off_energy.value) == 300.0
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        c.basis_precision = 'MEDIUM'
        assert issubclass(w[-1].category, UserWarning)
        assert "conflicts" in str(w[-1].message)
        assert c.param.cut_off_energy.value is None
        assert c.param.basis_precision.value.strip() == 'MEDIUM'

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        c.cut_off_energy = 200.0
        assert c.param.basis_precision.value is None
        assert issubclass(w[-1].category, UserWarning)
        assert 'option "cut_off_energy" conflicts' in str(w[-1].message)

    # test kpoint setup options
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # This block of tests is going to generate a lot of conflict warnings.
        # We already tested that those work, so just hide them from the output.

        c = Castep(kpts=[(0.0, 0.0, 0.0, 1.0),])
        assert c.cell.kpoint_list.value == '0.0 0.0 0.0 1.0'
        c.set_kpts(((0.0, 0.0, 0.0, 0.25), (0.25, 0.25, 0.3, 0.75)))
        assert c.cell.kpoint_list.value == '0.0 0.0 0.0 0.25\n0.25 0.25 0.3 0.75'
        c.set_kpts(c.cell.kpoint_list.value.split('\n'))
        assert c.cell.kpoint_list.value == '0.0 0.0 0.0 0.25\n0.25 0.25 0.3 0.75'
        c.set_kpts([3, 3, 2])
        assert c.cell.kpoint_mp_grid.value == '3 3 2'
        c.set_kpts(None)
        assert c.cell.kpoints_list.value is None
        assert c.cell.kpoint_list.value is None
        assert c.cell.kpoint_mp_grid.value is None
        c.set_kpts('2 2 3')
        assert c.cell.kpoint_mp_grid.value == '2 2 3'
        c.set_kpts({'even': True, 'gamma': True})
        assert c.cell.kpoint_mp_grid.value == '2 2 2'
        assert c.cell.kpoint_mp_offset.value == '0.25 0.25 0.25'
        c.set_kpts({'size': (2, 2, 4), 'even': False})
        assert c.cell.kpoint_mp_grid.value == '3 3 5'
        assert c.cell.kpoint_mp_offset.value == '0.0 0.0 0.0'
        atoms = ase_koopmans.build.bulk('Ag')
        atoms.calc = c
        c.set_kpts({'density': 10, 'gamma': False, 'even': None})
        assert c.cell.kpoint_mp_grid.value == '27 27 27'
        assert c.cell.kpoint_mp_offset.value == '0.018519 0.018519 0.018519'
        c.set_kpts({'spacing': (1 / (np.pi *10)), 'gamma': False, 'even': True})
        assert c.cell.kpoint_mp_grid.value == '28 28 28'
        assert c.cell.kpoint_mp_offset.value == '0.0 0.0 0.0'

    # test band structure setup
    from ase_koopmans.dft.kpoints import BandPath
    atoms = ase_koopmans.build.bulk('Ag')
    bp = BandPath(cell=atoms.cell,
                  path='GX',
                  special_points={'G': [0, 0, 0], 'X': [0.5, 0, 0.5]})
    bp = bp.interpolate(npoints=10)
    c = Castep(bandpath=bp)
    kpt_list = c.cell.bs_kpoint_list.value.split('\n')
    assert len(kpt_list) == 10
    assert list(map(float, kpt_list[0].split())) == [0., 0., 0.]
    assert list(map(float, kpt_list[-1].split())) == [0.5, 0.0, 0.5]
