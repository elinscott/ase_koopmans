def test_vasp2_import():
    """
    Test if we can find vasp2 using get_calculator()
    """

    from ase_koopmans.test.vasp import installed2 as installed
    from ase_koopmans.calculators.calculator import get_calculator_class

    assert installed()

    get_calculator_class('vasp2')
