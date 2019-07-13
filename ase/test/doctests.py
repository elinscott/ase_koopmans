import sys
from unittest import SkipTest


if sys.version_info < (3, 6):
    raise SkipTest('Test requires Python 3.6+, this is {}'
                   .format(sys.version_info))


import doctest
import importlib

module_names = """\
ase.atoms
ase.build.tools
ase.cell
ase.collections.collection
ase.dft.kpoints
ase.eos
ase.formula
ase.geometry.cell
ase.geometry.geometry
ase.io.ulm
ase.lattice
ase.spacegroup.findsym
ase.spacegroup.spacegroup
ase.spacegroup.xtal
ase.symbols
"""


for modname in module_names.splitlines():
    mod = importlib.import_module(modname)
    print(mod, doctest.testmod(mod, raise_on_error=True))
