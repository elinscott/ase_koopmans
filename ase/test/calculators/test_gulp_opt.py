def test_gulp_opt():
    import numpy as np
    from ase.calculators.gulp import GULP
    from ase.optimize import BFGS
    from ase.build import molecule, bulk
    from ase.constraints import ExpCellFilter

    # GULP optmization test
    atoms = molecule('H2O')
    atoms1 = atoms.copy()
    atoms1.calc = GULP(library='reaxff.lib')
    opt1 = BFGS(atoms1,trajectory='bfgs.traj')
    opt1.run(fmax=0.005)

    atoms2 = atoms.copy()
    calc2 = GULP(keywords='opti conp', library='reaxff.lib')
    opt2 = calc2.get_optimizer(atoms2)
    opt2.run()

    print(np.abs(opt1.atoms.positions - opt2.atoms.positions))
    assert np.abs(opt1.atoms.positions - opt2.atoms.positions).max() < 1e-5

    # GULP optimization test using stress
    atoms = bulk('Au', 'bcc', a=2.7, cubic=True)
    atoms1 = atoms.copy()
    atoms1.calc = GULP(keywords='conp gradient stress_out', library='reaxff_general.lib')
    atoms1f = ExpCellFilter(atoms1)
    opt1 = BFGS(atoms1f, trajectory='bfgs.traj')
    opt1.run(fmax=0.005)

    atoms2 = atoms.copy()
    calc2 = GULP(keywords='opti conp', library='reaxff_general.lib')
    opt2 = calc2.get_optimizer(atoms2)
    opt2.run()

    print(np.abs(opt1.atoms.positions - opt2.atoms.positions))
    assert np.abs(opt1.atoms.positions - opt2.atoms.positions).max() < 1e-5
