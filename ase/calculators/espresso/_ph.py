

""" Ph calculator

export ASE_PH_COMMAND="/path/to/ph PREFIX.phi"

"""


from ase import io
from ._espresso import EspressoParent


class EspressoPh(EspressoParent):
    """
    """
    ext_in = '.phi'
    ext_out = '.pho'
    implemented_properties = []
    command = 'ph.x -in PREFIX.phi > PREFIX.pho 2>&1'

    def __init__(self, *args, **kwargs):
        """
        All options for ph are copied verbatim to the input file
        """
        kwargs['label'] = 'ph'
        super().__init__(*args, **kwargs)
