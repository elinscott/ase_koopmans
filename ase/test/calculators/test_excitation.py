import numpy as np
import pytest
from ase.calculators.excitation import polarizability
from ase.calculators.h2morse import H2Morse
from ase.calculators.h2morse import H2MorseExcitedStates


def test_polarizabilty():
    """Test evaluation of polarizabily"""
    atoms = H2Morse()
    exl = H2MorseExcitedStates(atoms.get_calculator())

    

    alpha = polarizability(exl, range(2))
    assert alpha.shape == (2, )
    assert alpha.dtype == float
    alpha = polarizability(exl, 5 + 2j, tensor=True)
    assert alpha.shape == (3, 3)
    assert alpha.dtype == complex
    alpha = polarizability(exl, range(2), tensor=True)
    assert alpha.shape == (2, 3, 3)
    assert alpha.dtype == float


def main():
    test_polarizabilty()

if __name__ == '__main__':
    main()
