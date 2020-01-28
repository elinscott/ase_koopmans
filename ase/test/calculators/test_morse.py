from ase import Atoms
from ase.calculators.morse import MorsePotential

D0 = 5.
R0 = 3.


def test_gs_minimum_energy():
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, R0]])
    atoms.set_calculator(MorsePotential(epsilon=D0, rho=R0))
    assert atoms.get_potential_energy() == -D0

def test_gs_vibrations():
    # check ground state vibrations
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, R0]])
    atoms.set_calculator(MorsePotential(epsilon=D0, rho=R0))
    vib = Vibrations(atoms)
    vib.run()
    assert (vib.get_frequencies().real[-1] ==
            pytest.approx(ome[0], 1e-2))


test_minimum_energy()
