""" Calculator for koopmans_ham.x, a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kc_ham.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .khi and .kho

Run kc_ham.x jobs.
"""


import numpy as np
from ase.dft.kpoints import BandPath
from ._espresso import EspressoParent


class KoopmansHam(EspressoParent):
    ext_in = '.khi'
    ext_out = '.kho'
    implemented_properties = []

    # Default command does not use parallelism and assumes kc_ham.x is on the user's path
    command = 'kc_ham.x -in PREFIX.khi > PREFIX.kho'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_ham'
        super().__init__(*args, **kwargs)
