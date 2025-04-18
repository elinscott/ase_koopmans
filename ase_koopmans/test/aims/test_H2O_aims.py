from ase_koopmans import Atoms
from ase_koopmans.calculators.aims import Aims, AimsCube
from ase_koopmans.optimize import QuasiNewton


def test_H2O_aims():
    water = Atoms('HOH', [(1, 0, 0), (0, 0, 0), (0, 1, 0)])

    water_cube = AimsCube(points=(29, 29, 29),
                          plots=('total_density',
                                 'delta_density',
                                 'eigenstate 5',
                                 'eigenstate 6'))

    calc = Aims(xc='PBE',
                output=['dipole'],
                sc_accuracy_etot=1e-6,
                sc_accuracy_eev=1e-3,
                sc_accuracy_rho=1e-6,
                sc_accuracy_forces=1e-4,
                cubes=water_cube)

    water.calc = calc
    dynamics = QuasiNewton(water, trajectory='square_water.traj')
    dynamics.run(fmax=0.01)
