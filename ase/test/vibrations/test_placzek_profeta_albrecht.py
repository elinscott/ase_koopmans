"""
Test resonant Raman implementations
"""
import pytest

from ase.vibrations.placzek import Placzek, Profeta
from ase.vibrations.albrecht import Albrecht
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


def test_compare_placzek_implementation_intensities():
    """Intensities of different Placzek implementations
    should be similar"""
    atoms = H2Morse()
    name = 'placzek'
    pz = Placzek(atoms, H2MorseExcitedStates,
                 gsname=name, exname=name, txt=None)
    pz.run()
    om = 1
    pzi = pz.absolute_intensity(omega=om)[-1]

    # Profeta using frozenset
    pr = Profeta(atoms, H2MorseExcitedStates, approximation='Placzek',
                 gsname=name, exname=name, txt=None)
    pri = pr.absolute_intensity(omega=om)[-1]
    assert pzi == pytest.approx(pri, 1e-3)
    
    # Profeta using overlap
    name = 'profeta'
    pr = Profeta(atoms, H2MorseExcitedStates, approximation='Placzek',
                 gsname=name, exname=name,
                 overlap=lambda x, y: x.overlap(y),
                 txt=None)
    pr.run()
    pro = pr.absolute_intensity(omega=om)[-1]
    # print('pri, pro, pzi', pri, pro, pzi)
    assert pro == pytest.approx(pri, 1e-3)

def test_compare_placzek_albrecht_intensities():
    """Intensities of Placzek and Albrecht should be similar"""
    atoms = H2Morse()
    name = 'profeta'
    pr = Profeta(atoms, H2MorseExcitedStates, approximation='Placzek',
                 gsname=name, exname=name,
                 overlap=lambda x, y: x.overlap(y),
                 txt=None)
    pr.run()

    om = 10

    pr.approximation = 'p-p'
    pri = pr.absolute_intensity(omega=om)[-1]
    al = Albrecht(atoms, H2MorseExcitedStates,
                  gsname=name, exname=name,
                  approximation='Albrecht A', txt=None)
    ali = al.absolute_intensity(omega=om)[-1]
    print('pri, ali', pri, ali)
    
    """Albrecht and Placzek are approximately equal"""
    
    if 1:
        pri = pr.absolute_intensity(omega=om)[-1]
        al = Albrecht(atoms, H2MorseExcitedStates,
                      gsname=name, exname=name, overlap=True,
                      approximation='Albrecht', txt=None)
        ali = al.absolute_intensity(omega=om)[-1]
        print('pri, ali', pri, ali)
        
        """Albrecht A and P-P are approximately equal"""
        
        pr.approximation = 'p-p'
        pri = pr.absolute_intensity(omega=om)[-1]
        al.approximation = 'Albrecht A'
        ali = al.absolute_intensity(omega=om)[-1]
        print('pri, ali', pri, ali)
        assert pri == pytest.approx(ali, 1e-1)

        """Albrecht B+C and Profeta are approximately equal"""
        
        pr.approximation = 'Profeta'
        pri = pr.absolute_intensity(omega=om)[-1]
        al.approximation = 'Albrecht BC'
        ali = al.absolute_intensity(omega=om)[-1]
        print('pri, ali', pri, ali)
        ## assert pri == pytest.approx(ali, 1e-3)

def main():
    #test_compare_placzek_implementation_intensities()
    test_compare_placzek_albrecht_intensities()

if __name__ == '__main__':
    main()
