import numpy as np
import pytest

from ase.build import bulk
from ase.calculators.lj import LennardJones
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
from ase.optimize.precon.lbfgs import PreconLBFGS
from ase.constraints import UnitCellFilter, ExpCellFilter

spglib = pytest.importorskip('spglib')


def setup_cell():
    # setup an bcc Al cell
    at_init = bulk('Al', 'bcc', a=2 / np.sqrt(3), cubic=True)
    F = np.array([[1, 0.1, 0], [0.1, 1, 0], [0, 0, 1]])
    at_init.cell = at_init.cell @ F

    F = np.eye(3)
    for k in range(3):
        l = list(range(3))
        l.remove(k)
        (i, j) = l
        R = np.eye(3)
        theta = 0.1 * (k + 1)
        R[i, i] = np.cos(theta)
        R[j, j] = np.cos(theta)
        R[i, j] = np.sin(theta)
        R[j, i] = -np.sin(theta)
        F = np.dot(F, R)
    at_rot = at_init.copy()
    at_rot.set_cell(at_rot.cell @ F, True)
    return at_init, at_rot


def symmetrized_optimisation(at_init, filter):
    at = at_init.copy()
    calc = LennardJones()
    at.set_calculator(calc)

    at_cell = filter(at)
    dyn = PreconLBFGS(at_cell, precon=None)
    print("Energy", at.get_potential_energy(), at.get_volume())
    dyn.run(steps=300, fmax=0.001)
    print("n_steps", dyn.get_number_of_steps())
    print("Energy", at.get_potential_energy(), at.get_volume())
    print("Final forces\n", at.get_forces())
    print("Final stress\n", at.get_stress())

    print("initial symmetry at 1e-6")
    di = check_symmetry(at_init, 1.0e-6, verbose=True)
    print("final symmetry at 1e-6")
    df = check_symmetry(at, 1.0e-6, verbose=True)
    return di, df


@pytest.fixture(params=[UnitCellFilter, ExpCellFilter])
def filter(request):
    return request.param


# without symmetrization
# print("########### NO SYMMETRIZATION #############")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_no_symmetrization(filter):
    at_init, at_rot = setup_cell()
    at_unsym = at_init.copy()
    di, df = symmetrized_optimisation(at_unsym, filter)
    assert di["number"] == 63
    assert df["number"] <= 63


# print("######### ROTATED #########")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_no_sym_rotated(filter):
    at_init, at_rot = setup_cell()
    at_unsym_rot = at_rot.copy()
    di, df = symmetrized_optimisation(at_unsym_rot, filter)
    assert di["number"] == 63
    assert df["number"] <= 63


# symmetrization, adjust_positions, not cell
# print("########### SYMMETRIZATION, ADJUST POS BUT NO CELL #############")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_adj_pot(filter):
    at_init, at_rot = setup_cell()
    at_sym_2 = at_init.copy()
    at_sym_2.set_constraint(
        FixSymmetry(at_sym_2, adjust_positions=True, adjust_cell=False))
    di, df = symmetrized_optimisation(at_sym_2, filter)
    assert di["number"] == 63
    assert df["number"] >= 63


# print("######### ROTATED #########")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_rot_adj_pot(filter):
    at_init, at_rot = setup_cell()
    at_sym_2_rot = at_init.copy()
    at_sym_2_rot.set_constraint(
        FixSymmetry(at_sym_2_rot, adjust_positions=True, adjust_cell=False))
    di, df = symmetrized_optimisation(at_sym_2_rot, filter)
    assert di["number"] == 63
    assert df["number"] >= 63


# symmetrization, adjust_positions and cell
# print("########### SYMMETRIZATION, ADJUST POS AND CELL #############")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_adj_cell(filter):
    at_init, at_rot = setup_cell()
    at_sym_3 = at_init.copy()
    at_sym_3.set_constraint(
        FixSymmetry(at_sym_3, adjust_positions=True, adjust_cell=True))
    di, df = symmetrized_optimisation(at_sym_3, filter)
    assert di["number"] == 63
    assert df["number"] >= 63


# print("######### ROTATED #########")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_rot_adj_cell(filter):
    at_init, at_rot = setup_cell()
    at_sym_3_rot = at_init.copy()
    at_sym_3_rot.set_constraint(
        FixSymmetry(at_sym_3_rot, adjust_positions=True, adjust_cell=True))
    di, df = symmetrized_optimisation(at_sym_3_rot, filter)
    assert di["number"] == 63
    assert df["number"] >= 63
