""" Calculator for wann_to_kc.x, a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/wann_to_kc.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .w2ki and .w2ko

Run wann_to_kc.x jobs.
"""


from ._espresso import EspressoParent


class Wann2KC(EspressoParent):
    ext_in = '.w2ki'
    ext_out = '.w2ko'
    implemented_properties = []

    # Default command does not use parallelism and assumes wann2kc.x is on the user's path
    command = 'wann2kc.x -in PREFIX.w2ki > PREFIX.w2ko'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_screen'
        super().__init__(*args, **kwargs)
