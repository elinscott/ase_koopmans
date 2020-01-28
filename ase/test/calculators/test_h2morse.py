import pytest
from ase.vibrations import Vibrations
from ase.calculators.h2morse import H2Morse, Re, epsilon, Ee, ome
from ase.calculators.h2lj import H2LJExcitedStates


def test_gs_minimum():
    """Test ground state minimum distance and energy"""
    atoms = H2Morse()
    assert atoms.get_distance(0, 1) == Re[0]
    print(atoms.get_potential_energy(), -Ee[0])
    assert atoms.get_potential_energy() == -Ee[0]
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


test_gs_minimum()
if 0:
    test_excited_state()
