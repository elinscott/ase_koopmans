def test_fixcom():
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.optimize import BFGS
    from ase_koopmans.constraints import FixCom
    from ase_koopmans.build import molecule

    atoms = molecule('H2O')
    atoms.center(vacuum=4)
    atoms.calc = EMT()
    cold = atoms.get_center_of_mass()
    atoms.set_constraint(FixCom())

    opt = BFGS(atoms)
    opt.run(steps=5)

    cnew = atoms.get_center_of_mass()

    assert max(abs(cnew - cold)) < 1e-8
