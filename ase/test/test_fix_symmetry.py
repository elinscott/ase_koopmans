import numpy as np
import pytest

from ase.atoms import Atoms
from ase.calculators.lj import LennardJones
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
from ase.optimize.precon.lbfgs import PreconLBFGS
from ase.constraints import UnitCellFilter, ExpCellFilter

# setup an fcc Al cell
def setup_cell():
    a = 2.0/np.sqrt(3.0)
    at_prim = Atoms('Al2', positions=[[0,0,0],[a/2.0, a/2.0, a/2.0]],
                    cell=[[a,0,0],[0,a,0],[0,0,a]], pbc=[True, True, True])
    at_init = at_prim #* [2,2,2]
    F = np.array([[1,0.1,0],[0.1,1,0],[0,0,1]])
    at_init.set_cell(np.dot(at_init.get_cell(), F))

    F = np.eye(3)
    for k in range(3):
        l=list(range(3))
        l.remove(k)
        (i,j) = l
        R = np.eye(3)
        R[i,i] = np.cos(0.1*(k+1))
        R[j,j] = np.cos(0.1*(k+1))
        R[i,j] = np.sin(0.1*(k+1))
        R[j,i] = -np.sin(0.1*(k+1))
        F = np.dot(F, R)
    at_rot = at_init.copy()
    at_rot.set_cell(np.dot(at_rot.get_cell(),F), True)
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
#print("########### NO SYMMETRIZATION #############")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_no_symmetrization(filter):
    at_init, at_rot = setup_cell()
    at_unsym = at_init.copy()
    di, df = symmetrized_optimisation(at_unsym, filter)
    assert di["number"] == 63
    assert ((filter is ExpCellFilter and df["number"] == 63) or
            (filter is UnitCellFilter and df["number"] == 11))


#print("######### ROTATED #########")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_no_sym_rotated(filter):
    at_init, at_rot = setup_cell()
    at_unsym_rot = at_rot.copy()
    di, df = symmetrized_optimisation(at_unsym_rot, filter)
    assert di["number"] == 63
    assert ((filter is ExpCellFilter and df["number"] == 11) or
            (filter is UnitCellFilter and df["number"] == 2))
# symmetrization, adjust_positions, not cell
#print("########### SYMMETRIZATION, ADJUST POS BUT NO CELL #############")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_adj_pot(filter):
    at_init, at_rot = setup_cell()
    at_sym_2 = at_init.copy()
    at_sym_2.set_constraint(FixSymmetry(at_sym_2, adjust_positions=True, adjust_cell=False))
    di, df = symmetrized_optimisation(at_sym_2, filter)
    assert di["number"] == 63
    assert df["number"] == 139

#print("######### ROTATED #########")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_rot_adj_pot(filter):
    at_init, at_rot = setup_cell()
    at_sym_2_rot = at_init.copy()
    at_sym_2_rot.set_constraint(FixSymmetry(at_sym_2_rot, adjust_positions=True, adjust_cell=False))
    di, df = symmetrized_optimisation(at_sym_2_rot, filter)
    assert di["number"] == 63
    assert df["number"] == 139

# symmetrization, adjust_positions and cell
#print("########### SYMMETRIZATION, ADJUST POS AND CELL #############")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_adj_cell(filter):
    at_init, at_rot = setup_cell()
    at_sym_3 = at_init.copy()
    at_sym_3.set_constraint(FixSymmetry(at_sym_3, adjust_positions=True, adjust_cell=True))
    di, df = symmetrized_optimisation(at_sym_3, filter)
    assert di["number"] == 63
    assert df["number"] == 139

#print("######### ROTATED #########")
@pytest.mark.filterwarnings('ignore:ASE Atoms-like input is deprecated')
@pytest.mark.filterwarnings('ignore:Armijo linesearch failed')
def test_sym_rot_adj_cell(filter):
    at_init, at_rot = setup_cell()
    at_sym_3_rot = at_init.copy()
    at_sym_3_rot.set_constraint(FixSymmetry(at_sym_3_rot, adjust_positions=True, adjust_cell=True))
    di, df = symmetrized_optimisation(at_sym_3_rot, filter)
    assert di["number"] == 63
    assert df["number"] == 139
