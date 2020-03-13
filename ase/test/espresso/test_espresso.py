"""Check that QE calculation can run."""

from ase.build import bulk
from ase.calculators.espresso import Espresso


def verify(calc):
    assert calc.get_fermi_level() is not None
    assert calc.get_ibz_k_points() is not None
    assert calc.get_eigenvalues(spin=0, kpt=0) is not None
    assert calc.get_number_of_spins() is not None
    assert calc.get_k_point_weights() is not None


def test_main(espresso_factory):
    atoms = bulk('Si')
    calc = espresso_factory.calc()
    atoms.calc = calc
    atoms.get_potential_energy()
    verify(calc)


def test_smearing(espresso_factory):
    atoms = bulk('Cu')
    input_data = {'system':{'occupations': 'smearing',
                            'smearing': 'fermi-dirac',
                            'degauss': 0.02}}
    calc = espresso_factory.calc(input_data=input_data)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    verify(calc)
