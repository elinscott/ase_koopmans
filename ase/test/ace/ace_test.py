from ase import Atoms
from ase.calculators.acemolecule import ACE

label = "test"
mol = Atoms('H2',[(0, 0, 0),(0, 0, 0.7)])
basic = [dict(Cell= '5.0')]
ace = ACE(label=label, BasicInformation=basic)
mol.set_calculator(ace)
mol.get_forces()
#mol.get_potential_energy()
