def test_forcecurve(plt):
    from ase_koopmans.build import bulk
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.utils.forcecurve import force_curve
    from ase_koopmans.md import VelocityVerlet
    from ase_koopmans.units import fs
    from ase_koopmans.io import read

    atoms = bulk('Au', cubic=True) * (2, 1, 1)
    atoms.calc = EMT()
    atoms.rattle(stdev=0.05)

    md = VelocityVerlet(atoms, timestep=12.0 * fs, trajectory='tmp.traj')
    md.run(steps=10)
    images = read('tmp.traj', ':')
    force_curve(images)

    # import pylab as plt
    # plt.show()
