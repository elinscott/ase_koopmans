import numpy as np
import pytest
from ase.vibrations import Vibrations
from ase.calculators.h2morse import H2Morse, H2MorseState, Re, De, ome, Etrans
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
    assert (vib.get_frequencies().real[-1] ==
            pytest.approx(ome[0], 1e-2))


def test_gs_io_overlap():
    """Test ground state IO and 'wave function' overlap"""
    atoms0 = H2Morse()
    calc0 = atoms0.get_calculator()
    fname = 'calc0'
    calc0.write(fname)
    calc1 = H2MorseState(fname)
    for wf0, wf1 in zip(calc0.wfs, calc1.wfs):
        assert wf0 == pytest.approx(wf1, 1e-5)
    
    atoms1 = H2Morse()
    ov = calc0.overlap(calc0)
    # own overlap is the unity matrix
    assert np.eye(4) == pytest.approx(calc0.overlap(calc0), 1e-8)
    # self and other - test on unitarity
    ov = calc0.overlap(atoms1.get_calculator())
    assert np.eye(4) == pytest.approx(ov.dot(ov.T), 1e-8)
    
    
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


def test_excited_io():
    """Check writing and reading"""
    fname = 'exlist.dat'
    atoms = H2Morse()
    exl1 = H2MorseExcitedStates(atoms.get_calculator())
    exl1.write(fname)

    exl2 = H2MorseExcitedStates(fname)
    for ex1, ex2 in zip(exl1, exl2):
        assert ex1.energy == pytest.approx(ex2.energy, 1e-3)
        assert ex1.mur == pytest.approx(ex2.mur, 1e-5)
        assert ex1.muv == pytest.approx(ex2.muv, 1e-5)


def main():
    test_gs_io_overlap()

if __name__ == '__main__':
    main()
