from ase.calculators.test import FreeElectrons
from ase.geometry.bravais import all_variants
from ase.dft.band_structure import calculate_band_structure
from ase.utils import workdir
from ase import Atoms
import matplotlib.pyplot as plt

def test():
    ax = plt.gca()

    calc = FreeElectrons(nvalence=1)
    for i, lat in enumerate(all_variants()):
        if lat.ndim == 2:
            break
        xid = '{:02d}.{}'.format(i, lat.variant.name)
        atoms = Atoms(cell=lat.tocell())
        atoms.calc = calc
        path = lat.bandpath(density=80)
        path.write('path.{}.json'.format(xid))
        calc.set(kpts=path.kpts)

        atoms.calc = FreeElectrons(nvalence=1, kpts=path.kpts)

        bs = calculate_band_structure(atoms, path)
        bs.write('bs.{}.json'.format(xid))
        bs.plot(ax=ax, emin=-1, emax=22, filename='fig.{}.png'.format(xid))
        ax.clear()

with workdir('files', mkdir=True):
    test()
