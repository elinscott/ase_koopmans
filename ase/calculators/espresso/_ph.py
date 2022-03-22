## MARIJA TO UPDATE

""" wann2kcp calculator

export ASE_WANN2KCP_COMMAND="/path/to/wann2kcp PREFIX.wki"

"""


from ase import io
from ._espresso import EspressoParent


class EspressoPh(EspressoParent):
    """
    """
    ext_in = '.phi'
    ext_out = '.pho'
    implemented_properties = []
    command = 'wann2kcp.x -in PREFIX.wki > PREFIX.wko'

    def __init__(self, *args, **kwargs):
        """
        All options for wann2kcp are copied verbatim to the input file
        """
        kwargs['label'] = 'w2kcp'
        super().__init__(*args, **kwargs)
