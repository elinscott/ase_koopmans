""" Calculator for kcw.x (ham mode), a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kcw.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .khi and .kho

Run kcw.x jobs in the 'ham' mode.
"""


import numpy as np
from ase.dft.kpoints import BandPath
from ._espresso import EspressoParent


class KoopmansHam(EspressoParent):
    ext_in = '.khi'
    ext_out = '.kho'
    implemented_properties = []

    # Default command does not use parallelism and assumes kcw.x is on the user's path
    command = 'kcw.x -in PREFIX.khi > PREFIX.kho 2>&1'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_ham'
        super().__init__(*args, **kwargs)
