"""
Test resonant Raman implementations
"""
import pytest
import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.vibrations import Vibrations

from ase.vibrations.placzek import Placzek, Profeta
from ase.vibrations.albrecht import Albrecht
from ase.calculators.excitations import ExcitationList, Excitation
from ase.calculators.h2morse import H2Morse, H2MorseExcitedStates


def test_placzek_run():
    atoms = H2Morse()
    name = 'placzek'
    pz = Placzek(atoms, H2MorseExcitedStates,
                 gsname=name, exname=name, txt='-')
    pz.run()

def test_profeta_run():
    atoms = H2Morse()
    name = 'profeta'
    pr = Profeta(atoms, H2MorseExcitedStates,
                 gsname=name, exname=name, txt='-')
    pr.run()

def test_compare_intensities():
    atoms = H2Morse()
    pzname = 'placzek'
    pz = Placzek(atoms, H2MorseExcitedStates,
                 gsname=pzname, exname=pzname, txt=None)
    prname = 'profeta'
    prname = pzname
    pr = Profeta(atoms, H2MorseExcitedStates,
                 gsname=prname, exname=prname, txt=None)
    om = 1
    pzi = pz.absolute_intensity(omega=om)[-1]
    pri = pr.absolute_intensity(omega=om)[-1]
    pz.summary(om)
    print(pzi, pri)
    assert pzi == pytest.approx(pri, 1e-5)


def main():
    if 0:
        test_placzek_run()
        test_profeta_run()
    test_compare_intensities()

if __name__ == '__main__':
    main()
