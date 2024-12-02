def test_n2():
    from ase_koopmans import Atoms
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.optimize import QuasiNewton

    n2 = Atoms('N2', positions=[(0, 0, 0), (0, 0, 1.1)],
               calculator=EMT())
    QuasiNewton(n2).run(0.01)
    print(n2.get_distance(0, 1), n2.get_potential_energy())
