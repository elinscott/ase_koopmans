def test_vasp_cell():
    """

    Check the unit cell is handled correctly

    """

    from ase_koopmans.calculators.vasp import Vasp
    from ase_koopmans.build import molecule
    from ase_koopmans.test import must_raise

    # Molecules come with no unit cell

    atoms = molecule('CH4')
    calc = Vasp()

    with must_raise(RuntimeError):
        atoms.write('POSCAR')

    with must_raise(ValueError):
        atoms.calc = calc
        atoms.get_total_energy()
