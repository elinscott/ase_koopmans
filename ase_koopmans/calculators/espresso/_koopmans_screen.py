""" Calculator for kcw.x (screen mode), a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kcw.x -in PREFIX.ksi > PREFIX.kso"

N.B. the extensions must be .ksi and .kso

Run kcw.x jobs in the 'screen' mode.
"""


from ._espresso import EspressoParent


class KoopmansScreen(EspressoParent):
    ext_in = '.ksi'
    ext_out = '.kso'
    implemented_properties = []

    # Default command does not use parallelism and assumes kcw.x is on the user's path
    command = 'kcw.x -in PREFIX.ksi > PREFIX.kso 2>&1'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_screen'
        super().__init__(*args, **kwargs)
