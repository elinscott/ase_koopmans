def test_vasp2_cell():
    """

    Check the unit cell is handled correctly

    """

    from ase_koopmans.test.vasp import installed2 as installed
    from ase_koopmans.calculators.vasp import Vasp2 as Vasp
    from ase_koopmans.build import molecule
    from ase_koopmans.test import must_raise
    assert installed()


    # Molecules come with no unit cell

    atoms = molecule('CH4')
    calc = Vasp()

    with must_raise(ValueError):
        atoms.calc = calc
        atoms.get_total_energy()
