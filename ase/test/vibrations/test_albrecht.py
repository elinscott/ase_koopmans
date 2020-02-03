import numpy as np
import pytest
from ase.calculators.h2morse import H2Morse, H2MorseExcitedStates
from ase.vibrations.placzek import Placzek, Profeta
from ase.vibrations.albrecht import Albrecht


def test_overlap():
    name = 'rrmorse'
    atoms = H2Morse()
    om = 1

    ao = Albrecht(atoms, H2MorseExcitedStates,
                  gsname=name, exname=name, 
                  overlap=lambda x, y: x.overlap(y),
                  approximation='Albrecht A', txt=None)
    ao.run()

    """One state only"""
    
    ao = Albrecht(atoms, H2MorseExcitedStates, exkwargs={'nstates':1},
                  gsname=name, exname=name, overlap=True,
                  approximation='Albrecht A', txt=None)
    ao.run()
    ao = Albrecht(atoms, H2MorseExcitedStates, exkwargs={'nstates':1},
                  gsname=name, exname=name, overlap=True,
                  approximation='Albrecht A', txt=None)
    ao.run()
    aoi = ao.absolute_intensity(omega=om)[-1]
    al = Albrecht(atoms, H2MorseExcitedStates, exkwargs={'nstates':1},
                  gsname=name, exname=name, 
                  approximation='Albrecht A', txt=None)
    ali = al.absolute_intensity(omega=om)[-1]
    print('w/o and with overlap', ali, aoi)
    assert ali == pytest.approx(aoi, 1e-9)

    """Include degenerate states"""
    
    ao = Albrecht(atoms, H2MorseExcitedStates, exkwargs={'nstates':1},
                  gsname=name, exname=name, overlap=True,
                  approximation='Albrecht A', txt=None)
    ao.run()
    aoi = ao.absolute_intensity(omega=om)[-1]
    al = Albrecht(atoms, H2MorseExcitedStates, exkwargs={'nstates':1},
                  gsname=name, exname=name, 
                  approximation='Albrecht A', txt=None)
    ali = al.absolute_intensity(omega=om)[-1]
    print('w/o and with overlap', ali, aoi)
    assert ali == pytest.approx(aoi, 1e-9)


def test_displacements():
    name = 'albrecht'
    atoms = H2Morse()
    al = Albrecht(atoms, H2MorseExcitedStates,
                  gsname=name, exname=name, 
                  approximation='Albrecht A', txt=None)
    al.run()
    al.summary()

    pr = Profeta(atoms, H2MorseExcitedStates, approximation='P-P',
                 gsname=name, exname=name,
                 txt=None)
    pr.summary()

    
def main():
    #test_compare_placzek_implementation_intensities()
    test_overlap()

if __name__ == '__main__':
    main()
