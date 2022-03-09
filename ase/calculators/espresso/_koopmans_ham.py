""" Calculator for kcw.x (ham mode), a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kcw.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .khi and .kho

Run kcw.x jobs in the 'ham' mode.
"""


import numpy as np
from ase.dft.kpoints import BandPath
from ._espresso import EspressoParent, EspressoWithBandstructure


class KoopmansHam(EspressoWithBandstructure, EspressoParent):
    ext_in = '.khi'
    ext_out = '.kho'
    implemented_properties = []

    # Default command does not use parallelism and assumes kcw.x is on the user's path
    command = 'kcw.x -in PREFIX.khi > PREFIX.kho'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_ham'
        super().__init__(*args, **kwargs)

    def read_results(self):
        super().read_results()
        if isinstance(self.parameters.kpts, BandPath) and len(self.parameters.kpts.kpts) > 1:
            # Add the bandstructure to the results
            self.band_structure(vbm_to_zero=True)

    @property
    def vbm_index(self):
        return self.parameters.num_wann_occ - 1

    def eigenvalues_from_results(self):
        assert 'eigenvalues' in self.results, 'Please call {0}.calculate() prior to calling {0}.band_structure'.format(
            self.__class__.__name__)

        return np.array([self.results['eigenvalues']])
