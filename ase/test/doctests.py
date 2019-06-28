import doctest
import sys

from ase import atoms, cell
from ase.collections import collection
from ase import lattice
from ase.dft import kpoints
from ase.spacegroup import spacegroup, findsym, xtal
from ase.geometry import geometry, cell as geometry_cell
from ase.build import tools
from ase.io import ulm
import ase.eos as eos
import ase.formula as formula

modules = [xtal, spacegroup, geometry_cell, findsym, atoms, eos,
           geometry, tools, collection]

if sys.version_info >= (3, 6):
    # dicts must be ordered
    modules += [ulm, formula, kpoints, cell, lattice]

for mod in modules:
    print(mod, doctest.testmod(mod, raise_on_error=True))
