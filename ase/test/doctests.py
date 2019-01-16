import doctest

from ase import atoms
from ase.collections import collection
from ase.spacegroup import spacegroup, findsym, xtal
from ase.geometry import geometry, cell
from ase.build import tools
from ase.io import ulm
import ase.eos as eos

modules = [xtal, spacegroup, cell, findsym, ulm, atoms, eos,
           geometry, tools, collection]

for mod in modules:
    print(mod, doctest.testmod(mod, raise_on_error=True))
