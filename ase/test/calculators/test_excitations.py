import pytest
from ase.calculators.excitations import polarizability
from ase.calculators.h2morse import H2Morse
from ase.calculators.h2morse import H2MorseExcitedStates


def test_polarizabilty():
    """Test evaluation of polarizabily"""
    atoms = H2Morse()
    exl = H2MorseExcitedStates(atoms.get_calculator())

    alphaf = polarizability(exl, range(2))
    assert alphaf.shape == (2, )
    assert alphaf.dtype == float
    alphat = polarizability(exl, 5 + 2j, tensor=True)
    assert alphat.shape == (3, 3)
    assert alphat.dtype == complex
    alphat = polarizability(exl, range(2), tensor=True)
    assert alphat.shape == (2, 3, 3)
    assert alphat.dtype == float

    # check tensor
    for af, at in zip(alphaf, alphat):
        assert at.diagonal().sum() / 3 == pytest.approx(af, 1.e-8)


def main():
    test_polarizabilty()

if __name__ == '__main__':
    main()
