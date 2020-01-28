import pytest
from ase.vibrations import Vibrations
from ase.calculators.h2morse import H2Morse, Re, De, ome, Etrans
from ase.calculators.h2morse import H2MorseExcitedStates


def test_gs_minimum():
    """Test ground state minimum distance, energy and 
    vibrational frequency"""
    atoms = H2Morse()
    assert atoms.get_distance(0, 1) == Re[0]
    assert atoms.get_potential_energy() == -De[0]
    # check ground state vibrations
    vib = Vibrations(atoms)
    vib.run()
    #vib.summary()
    assert (vib.get_frequencies().real[-1] ==
            pytest.approx(ome[0], 1e-2))


def test_excited_state():
    """Test excited state transition energies"""
    gsatoms = H2Morse()
    Egs0 = gsatoms.get_potential_energy()
    for i in range(1, 4):
        exatoms = H2Morse()
        exatoms[1].position[2] = Re[i]
        Egs = exatoms.get_potential_energy()
        exl = H2MorseExcitedStates(exatoms.get_calculator())
        assert (exl[i - 1].energy ==
                pytest.approx(Etrans[i] - Egs + Egs0, 1e-8))

