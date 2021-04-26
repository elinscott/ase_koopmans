""" Calculator for koopmans_screen.x, a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kc_screen.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .ksi and .kso

Run kc_screen.x jobs.
"""


from ._espresso import EspressoParent


class KoopmansScreen(EspressoParent):
    ext_in = 'ksi'
    ext_out = 'kso'
    implemented_properties = []

    # Default command does not use parallelism and assumes kc_screen.x is on the user's path
    command = 'kc_screen.x -in PREFIX.ksi > PREFIX.kso'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_screen'
        super().__init__(*args, **kwargs)
