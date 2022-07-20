from ase import io
from ._espresso import EspressoParent


class Projwfc(EspressoParent):
    """
    """
    ext_in = '.pri'
    ext_out = '.pro'
    implemented_properties = []
    command = 'projwfc.x -in PREFIX.pri > PREFIX.pro 2>&1'

    def __init__(self, *args, **kwargs):
        """
        All options for projwfc are copied verbatim to the input file
        """
        kwargs['label'] = 'projwfc'
        super().__init__(*args, **kwargs)
