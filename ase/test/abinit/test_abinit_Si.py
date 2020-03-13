def test_abinit_Si(abinit_factory):
    from ase.build import bulk
    from ase.units import Ry
    from ase.calculators.abinit import Abinit

    atoms = bulk('Si')

    calc = abinit_factory.calc(
        label='Si',
        nbands=8,
        ecut=10 * Ry,
        kpts=[4, 4, 4],
    )

    calc.set(toldfe=1.0e-2)
    atoms.set_calculator(calc)
    e = atoms.get_potential_energy()
    print(e)
