""" pw2wannier calculator

export ASE_PW2WANNIER_COMMAND="/path/to/pw2wannier PREFIX.p2wi"
                                                       
"""


from ase import io
from ase.calculators.calculator import FileIOCalculator


class PW2Wannier(FileIOCalculator):
    """
    """
    implemented_properties = []
    command = 'pw2wannier90.x PREFIX.p2wi > PREFIX.p2wo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='pw2wann', atoms=None, **kwargs):
        """
        All options for pw2wannier are copied verbatim to the input file

        """
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.calc = None

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        # Create the appropriate directory
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        # Write the input file
        io.write(self.label + '.p2wi', atoms)

    def read_results(self):
        output = io.read(self.label + '.p2wo')
        self.calc = output.calc
        self.results = output.calc.results
