import pytest
import numpy as np

from ase.build import molecule
from ase.build.attach import attach, attach_randomly


def test_attach():
    """Attach two molecules and check that their minimal distance
    is as required"""
    m1 = molecule('C6H6')
    m2 = molecule('NH3')

    distance = 2.
    m12 = attach(m1, m2, distance)
    dmin = np.linalg.norm(m12[15].position - m12[8].position)
    assert dmin == pytest.approx(distance, 1e-8)


def test_attach_randomly():
    m1 = molecule('C6H6')
    m2 = molecule('CF4')

    np.random.seed(42)
    
    distance = 2.5
    pos2_ac = np.zeros((5, 3))
    N = 25
    for i in range(N):
        atoms = attach_randomly(m1, m2, distance)
        pos2_ac += atoms.get_positions()[12:, :]
    # the average position should be "zero" approximately
    assert (np.abs(pos2_ac / N) <= 1).all()


if __name__ == '__main__':
    test_attach_randomly()
