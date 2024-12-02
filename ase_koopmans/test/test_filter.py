def test_filter():
    """Test that the filter and trajectories are playing well together."""

    from ase_koopmans.build import molecule
    from ase_koopmans.constraints import Filter
    from ase_koopmans.optimize import QuasiNewton
    from ase_koopmans.calculators.emt import EMT

    atoms = molecule('CO2')
    atoms.calc = EMT()
    filter = Filter(atoms, indices=[1, 2])

    opt = QuasiNewton(filter, trajectory='filter-test.traj', logfile='filter-test.log')
    opt.run()
