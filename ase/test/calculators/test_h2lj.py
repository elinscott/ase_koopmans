import pytest
from ase.calculators.h2lj import H2LJ, Re, epsilon, Ee
from ase.calculators.h2lj import H2LJExcitedStates

def test_gs_minimum():
    """Test ground state minimum distance and energy"""
    atoms = H2LJ()
    assert atoms.get_distance(0, 1) == Re[0]
    assert atoms.get_potential_energy() == -epsilon[0]

def test_ex_state_energy():
    """Test excited state transition energies"""
    atoms = H2LJ()
    exlist = H2LJExcitedStates(atoms.get_calculator())
    for i, ex in enumerate(exlist):
        ## print(i, ex)
        assert ex.energy == pytest.approx(Ee[i + 1], 1e-8)

test_ex_state_energy()
