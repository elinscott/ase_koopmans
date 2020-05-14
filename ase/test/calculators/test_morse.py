from ase import Atoms
from ase.vibrations import Vibrations
from ase.calculators.morse import MorsePotential

De = 5.
Re = 3.
rho0 = 2.


def test_gs_minimum_energy():
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, Re]])
    atoms.calc = MorsePotential(epsilon=De, r0=Re)
    assert atoms.get_potential_energy() == -De


def test_gs_vibrations():
    # check ground state vibrations
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, Re]])
    atoms.calc = MorsePotential(epsilon=De, r0=Re, rho0=rho0)
    vib = Vibrations(atoms)
    vib.run()
