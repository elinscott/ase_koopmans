from ase import Atoms
from ase.calculators.lj import LennardJones

atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 2 ** (1.0 / 6.0)]])


def test_minimum_energy():
    # See https://en.wikipedia.org/wiki/Lennard-Jones_potential
    atoms.set_calculator(LennardJones(rc=1.0e5))
    assert atoms.get_potential_energy() == -1.0


def test_system_changes():
    atoms.set_calculator(LennardJones(rc=1.0e5))
    atoms.calc.calculate(atoms, system_changes=["positions"])
