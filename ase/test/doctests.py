import doctest
import sys

from ase import atoms
from ase.collections import collection
from ase.spacegroup import spacegroup, findsym, xtal
from ase.geometry import geometry, cell
from ase.build import tools
from ase.io import ulm
import ase.eos as eos
import ase.formula as formula

modules = [xtal, spacegroup, cell, findsym, atoms, eos,
           geometry, tools, collection]

if sys.version_info >= (3, 6):
    # dicts must be ordered
    modules.append(ulm)
    modules.append(formula)

for mod in modules:
    print(mod, doctest.testmod(mod, raise_on_error=True))
