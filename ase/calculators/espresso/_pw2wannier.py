""" pw2wannier calculator

export ASE_PW2WANNIER_COMMAND="/path/to/pw2wannier PREFIX.p2wi"

"""


from ase import io
from ._espresso import EspressoParent


class PW2Wannier(EspressoParent):
    """
    """
    ext_in = '.p2wi'
    ext_out = '.p2wo'
    implemented_properties = []
    command = 'pw2wannier90.x -in PREFIX.p2wi > PREFIX.p2wo 2>&1'

    def __init__(self, *args, **kwargs):
        """
        All options for pw2wannier are copied verbatim to the input file
        """
        kwargs['label'] = 'pw2wann'
        super().__init__(*args, **kwargs)
