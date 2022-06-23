""" Calculator for kcw.x (wann2kcw mode), a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/kcw.x -in PREFIX.w2ki > PREFIX.w2ko"

N.B. the extensions must be .w2ki and .w2ko

Run kcw.x jobs in the 'wann2kcw' mode.
"""


from ._espresso import EspressoParent


class Wann2KC(EspressoParent):
    ext_in = '.w2ki'
    ext_out = '.w2ko'
    implemented_properties = []

    # Default command does not use parallelism and assumes kcw.x is on the user's path
    command = 'kcw.x -in PREFIX.w2ki > PREFIX.w2ko 2>&1'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'wann2kc'
        super().__init__(*args, **kwargs)
