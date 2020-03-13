"""Check that QE calculation can run."""

from ase.build import bulk
from ase.calculators.espresso import Espresso


def test_main(espresso_factory):
    silicon = bulk('Si')
    calc = espresso_factory.calc()
    silicon.calc = calc
    print(silicon.calc.parameters)
    silicon.get_potential_energy()

    assert calc.get_fermi_level() is not None
    assert calc.get_ibz_k_points() is not None
    assert calc.get_eigenvalues(spin=0, kpt=0) is not None
    assert calc.get_number_of_spins() is not None
    assert calc.get_k_point_weights() is not None
