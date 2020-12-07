""" wannier90 calculator

export ASE_WANNIER90_COMMAND="/path/to/wannier90 PREFIX.win"

"""

from ase import io
from ase.calculators.calculator import FileIOCalculator

class Wannier90(FileIOCalculator):
    """
    """
    implemented_properties = []
    command = 'wannier90 PREFIX.pw2wan'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='wannier', atoms=None, **kwargs):
        """
        All options for wannier90 are copied verbatim to the input file

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
        io.write(self.label + '.win', atoms)

    def read_results(self):
        output = io.read(self.label + '.wout')
        self.calc = output.calc
        self.results = output.calc.results
