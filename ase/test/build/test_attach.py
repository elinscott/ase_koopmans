import pytest
import numpy as np

from ase.build import molecule
from ase.build.attach import attach


def test_attach():
    """Attach two molecules and check that their minimal distance
    is as required"""
    m1 = molecule('C6H6')
    m2 = molecule('NH3')

    distance = 2.
    m12 = attach(m1, m2, distance)
    dmin = np.linalg.norm(m12[15].position - m12[8].position)
    assert dmin == pytest.approx(distance, 1e-8)
    
    m12.edit()


if __name__ == '__main__':
    test_attach()
