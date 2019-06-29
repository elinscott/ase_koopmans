from ase.build import bulk
from ase.calculators.siesta import Siesta
from ase.utils import workdir

atoms = bulk('Si')
with workdir('files', mkdir=True):
    path = atoms.cell.bandpath('GXWK', density=10)
    calc = Siesta(kpts=[2, 2, 2], bandpath=path)
    atoms.calc = calc
    atoms.get_potential_energy()
    bs = atoms.calc.band_structure()
    print(bs)
    bs.write('bs.json')
