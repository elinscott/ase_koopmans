# creates: cu.png
from ase_koopmans.build import bulk
from ase_koopmans.calculators.test import FreeElectrons

a = bulk('Cu')
a.calc = FreeElectrons(nvalence=1,
                       kpts={'path': 'GXWLGK', 'npoints': 200})
a.get_potential_energy()
bs = a.calc.band_structure()
bs.plot(emax=10, filename='cu.png')
