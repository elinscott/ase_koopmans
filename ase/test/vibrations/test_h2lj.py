from ase.calculators.lj import LennardJones
from ase.calculators.h2_ljexcitation import H2LJ, Re, epsilon

def test_lj():
    atoms = H2LJ()
    calc = atoms.get_calculator()
    calc.parameters.epsilon = 1
    calc.parameters.sigma = 1
    atoms[1].position[2] = 2**(1. / 6.)
    print(atoms.get_potential_energy())
    assert atoms.get_potential_energy() == -1.0
    

def test_gs_minimum():
    atoms = H2LJ()
    print(atoms.get_calculator().parameters, Re[0], atoms.get_calculator().parameters['sigma'] * 2**(1/6))
    assert atoms.get_distance(0, 1) == Re[0]
    assert atoms.get_potential_energy() == -epsilon[0]


