def test_turbomole_H2():
    from ase_koopmans import Atoms
    from ase_koopmans.calculators.turbomole import Turbomole

    atoms = Atoms('H2', positions=[(0, 0, 0), (0, 0, 1.1)])

    # Write all commands for the define command in a string
    define_str = '\n\na coord\n*\nno\nb all sto-3g hondo\n*\neht\n\n\n\n*'

    atoms.calc = Turbomole(define_str=define_str)

    # Run turbomole
    atoms.get_potential_energy()
