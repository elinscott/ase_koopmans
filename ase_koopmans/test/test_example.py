def test_example():
    from ase_koopmans import Atoms
    from ase_koopmans.constraints import FixAtoms
    from ase_koopmans.io import Trajectory
    from ase_koopmans.optimize import QuasiNewton
    from ase_koopmans.calculators.morse import MorsePotential

    atoms = Atoms('H7',
                  positions=[(0, 0, 0),
                             (1, 0, 0),
                             (0, 1, 0),
                             (1, 1, 0),
                             (0, 2, 0),
                             (1, 2, 0),
                             (0.5, 0.5, 1)],
                  constraint=[FixAtoms(range(6))],
                  calculator=MorsePotential())

    traj = Trajectory('H.traj', 'w', atoms)
    dyn = QuasiNewton(atoms, maxstep=0.2)
    dyn.attach(traj.write)
    dyn.run(fmax=0.01, steps=100)

    print(atoms)
    del atoms[-1]
    print(atoms)
    del atoms[5]
    print(atoms)
    assert len(atoms.constraints[0].index) == 5
