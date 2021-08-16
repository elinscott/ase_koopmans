""" Calculator for koopmans_ham.x, a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kc_ham.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .khi and .kho

Run kc_ham.x jobs.
"""


from ._espresso import EspressoParent, EspressoWithBandstructure


class KoopmansHam(EspressoParent, EspressoWithBandstructure):
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
        self.band_structure(vbm_to_zero=True)

    @property
    def vbm_index(self):
        return self.parameters['input_data']['wannier']['num_wann_occ'] - 1

    def eigenvalues_from_results(self):
        assert 'eigenvalues' in self.results, 'Please call {0}.calculate() prior to calling {0}.band_structure'.format(
            self.__class__.__name__)

        return np.array([self.results['eigenvalues']])
