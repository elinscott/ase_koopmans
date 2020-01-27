from ase import Atoms
from ase.calculators.lj import LennardJones


def test_minimum_energy():
    # See https://en.wikipedia.org/wiki/Lennard-Jones_potential
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 2**(1. / 6.)]])
    atoms.set_calculator(LennardJones(rc=1.e5))
    assert atoms.get_potential_energy() == -1.0
