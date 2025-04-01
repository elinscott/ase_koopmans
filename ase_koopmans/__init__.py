# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""Atomic Simulation Environment."""

import sys

import numpy as np
from packaging.version import Version

if sys.version_info[0] == 2:
    raise ImportError('ASE requires Python3. This is Python2.')


if Version(np.__version__) < '1.9':
    raise ImportError(
        'ASE needs NumPy-1.9.0 or later. You have: %s' % np.__version__)


__all__ = ['Atoms', 'Atom']

__version__ = '0.1.8'

# import ase_koopmans.parallel early to avoid circular import problems when
# ase_koopmans.parallel does "from gpaw.mpi import world":
import ase_koopmans.parallel  # noqa
from ase_koopmans.atom import Atom
from ase_koopmans.atoms import Atoms

ase_koopmans.parallel  # silence pyflakes
