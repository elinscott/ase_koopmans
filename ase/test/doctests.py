import doctest
import importlib
import shutil
import sys
from unittest import SkipTest
from distutils.version import LooseVersion

import numpy as np

if sys.version_info < (3, 6):
    raise SkipTest('Test requires Python 3.6+, this is {}'
                   .format(sys.version_info))


# Older numpies format arrays differently.
# We use the printoptions contextmanager from numpy 1.15:
# https://docs.scipy.org/doc/numpy/release.html#id45
if LooseVersion(np.__version__) < '1.15':
    raise SkipTest('need numpy >= 1.15')


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
ase.phasediagram
ase.spacegroup.findsym
ase.spacegroup.spacegroup
ase.spacegroup.xtal
ase.symbols
"""


for modname in module_names.splitlines():
    if modname == 'ase.spacegroup.findsym' and not shutil.which('findsym'):
        print('Skipping {} because we do not have findsym'.format(modname))
        continue

    mod = importlib.import_module(modname)
    with np.printoptions(legacy='1.13'):
        print(mod, doctest.testmod(mod, raise_on_error=True))
