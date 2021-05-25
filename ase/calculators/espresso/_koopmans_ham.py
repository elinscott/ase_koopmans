""" Calculator for koopmans_ham.x, a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kc_ham.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .khi and .kho

Run kc_ham.x jobs.
"""


import numpy as np
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure
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

    def read_results(self):
        super().read_results()
        # Add the bandstructure to the results
        self.band_structure()

    def band_structure(self):
        # Construct bandstructure here (rather than within self.calculate()) because we have access to the band path
        assert 'eigenvalues' in self.results, 'Please call {0}.calculate() prior to calling {0}.band_structure'.format(
            self.__class__.__name__)

        eigenvalues_np = np.array([self.results['eigenvalues']])

        # Shift so that the VBM = 0.0
        try:
            num_occ = self.parameters['input_data']['wannier']['num_wann_occ']
            if num_occ > 0:
                reference = np.max(eigenvalues_np[:, :, num_occ - 1])
                eigenvalues_np -= reference
        except KeyError():
            pass

        self.results['band structure'] = BandStructure(self.parameters['kpts'], eigenvalues_np)
        return self.results['band structure']
