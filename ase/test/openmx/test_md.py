def test_md():
    import numpy as np
    import unittest

    from ase.units import Ry, Ha
    from ase.calculators.openmx import OpenMX
    from ase.io.trajectory import Trajectory
    from ase.optimize import QuasiNewton
    from ase.constraints import UnitCellFilter
    from ase.calculators.calculator import PropertyNotImplementedError
    from ase import Atoms
    """ Only OpenMX 3.8 or higher version pass this test"""

    bud = Atoms('CH4', np.array([
            [0.000000, 0.000000, 0.100000],
            [0.682793, 0.682793, 0.682793],
            [-0.682793, -0.682793, 0.68279],
            [-0.682793, 0.682793, -0.682793],
            [0.682793, -0.682793, -0.682793]]),
            cell=[10, 10, 10])

    calc = OpenMX(
        label='ch4',
        xc='GGA',
        energy_cutoff=300 * Ry,
        convergence=1e-4 * Ha,
        # Use 'C_PBE19' and 'H_PBE19' for version 3.9
        definition_of_atomic_species=[['C', 'C5.0-s1p1', 'C_PBE13'],
                                      ['H', 'H5.0-s1', 'H_PBE13']],
        kpts=(4, 4, 4),
        eigensolver='Band'
        )

    bud.calc = calc

    try:
        bud.get_stress()
    except PropertyNotImplementedError as err:
        raise unittest.SkipTest(err)

    traj = Trajectory('example.traj', 'w', bud)
    ucf = UnitCellFilter(bud, mask=[True, True, False, False, False, False])
    dyn = QuasiNewton(ucf)
    dyn.attach(traj.write)
    dyn.run(fmax=0.1)
    bud.get_potential_energy()
    # XXX maybe assert something?

    traj.close()
