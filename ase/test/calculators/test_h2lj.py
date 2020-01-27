import pytest
from ase.vibrations import Vibrations
from ase.calculators.h2lj import H2LJ, Re, epsilon, Ee, ome
from ase.calculators.h2lj import H2LJExcitedStates


def test_gs_minimum():
    """Test ground state minimum distance and energy"""
    atoms = H2LJ()
    assert atoms.get_distance(0, 1) == Re[0]
    assert atoms.get_potential_energy() == -epsilon[0]
    # check ground state vibrations
    vib = Vibrations(atoms)
    vib.run()
    assert (vib.get_frequencies().real[-1] ==
            pytest.approx(ome[0], 1e-2))


def test_excited_state():
    """Test excited state transition energies"""
    atoms = H2LJ()
    exlist = H2LJExcitedStates(atoms.get_calculator())
    for i, ex in enumerate(exlist):
        #print(i, ex)
        assert ex.energy == pytest.approx(Ee[i + 1], 1e-8)

if 0:
    test_gs_minimum()
    test_excited_state()
