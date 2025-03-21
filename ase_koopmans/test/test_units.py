"""This test cross-checks our implementation of CODATA against the
implementation that SciPy brings with it.
"""


def test_units():
    import numpy as np
    from ase_koopmans.units import CODATA
    import scipy.constants.codata

    name_map = {'_c': 'speed of light in vacuum',
                '_mu0': 'mag. const.',
                '_Grav': 'Newtonian constant of gravitation',
                '_hplanck': 'Planck constant',
                '_e': 'elementary charge',
                '_me': 'electron mass',
                '_mp': 'proton mass',
                '_Nav': 'Avogadro constant',
                '_k': 'Boltzmann constant',
                '_amu': 'atomic mass unit-kilogram relationship'}

    for version in sorted(CODATA.keys()):
        print('Checking CODATA version "{0}"'.format(version))

        try:
            scipy_CODATA = getattr(scipy.constants.codata,
                                   '_physical_constants_{0}'.format(version))
        except AttributeError:
            print('\tNot available through scipy, skipping')
            continue

        for unit, scipyname in name_map.items():
            ase_koopmansval = CODATA[version][unit]
            try:
                scipyval = scipy_CODATA[name_map[unit]][0]
                msg = 'Unit "{0}" : '.format(name_map[unit])
                ok = True
                if np.isclose(ase_koopmansval, scipyval):
                    msg += '[OK]'
                else:
                    msg += '[FALSE]'
                    ok = False
                print('\t' + msg)
                if not ok:
                    raise AssertionError

            except KeyError:
                # 2002 in scipy contains too little data
                continue


def test_create_units():
    """Check that units are created and allow attribute access."""

    import ase_koopmans.units

    print('Checking create_units and attribute access')
    # just use current CODATA version
    new_units = ase_koopmans.units.create_units(ase_koopmans.units.__codata_version__)
    assert new_units.eV == new_units['eV'] == ase_koopmans.units.eV
    for unit_name in new_units.keys():
        assert getattr(new_units, unit_name) == getattr(ase_koopmans.units, unit_name)
        assert new_units[unit_name] == getattr(ase_koopmans.units, unit_name)
