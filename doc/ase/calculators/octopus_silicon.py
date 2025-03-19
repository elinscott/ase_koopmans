from ase_koopmans.build import bulk
from ase_koopmans.calculators.octopus import Octopus

system = bulk('Si', orthorhombic=True)

calc = Octopus(label='silicon',
               Spacing=0.25,
               KPointsGrid=[[4, 4, 4]],
               KPointsUseSymmetries=True,
               Output='dos + density + potential',
               OutputFormat='xcrysden',
               DosGamma=0.1)

system.calc = calc
system.get_potential_energy()
