from ase_koopmans.build import bulk
from ase_koopmans.calculators.dftd3 import DFTD3

diamond = bulk('C')
d3 = DFTD3()
diamond.calc = d3
diamond.get_potential_energy()
