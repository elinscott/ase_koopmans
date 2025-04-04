def test_fix_bond_length_mic():
    import ase_koopmans
    from ase_koopmans.calculators.lj import LennardJones
    from ase_koopmans.constraints import FixBondLength
    from ase_koopmans.optimize import FIRE

    for wrap in [False, True]:
        a = ase_koopmans.Atoms('CCC',
                      positions=[[1, 0, 5],
                                 [0, 1, 5],
                                 [-1, 0.5, 5]],
                      cell=[10, 10, 10],
                      pbc=True)

        if wrap:
            a.set_scaled_positions(a.get_scaled_positions() % 1.0)
        a.calc = LennardJones()
        a.set_constraint(FixBondLength(0, 2))

        d1 = a.get_distance(0, 2, mic=True)

        FIRE(a, logfile=None).run(fmax=0.01)
        e = a.get_potential_energy()
        d2 = a.get_distance(0, 2, mic=True)
        assert abs(e - -2.034988) < 1e-6
        assert abs(d1 - d2) < 1e-6
