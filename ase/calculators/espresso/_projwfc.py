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
    
    def read_xml(self):
        """
        Read the atomic_proj.xml file and store its contents in self.results
        """
        xml_file = self.parameters['outdir'] / (self.parameters['prefix'] + '.save') / 'atomic_proj.xml'
        with open(xml_file, 'r') as fd:
            xml_dict = io.espresso.read_projwfc_xml(fd)

        self.results.update(**xml_dict)

        return 