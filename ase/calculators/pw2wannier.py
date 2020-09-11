""" pw2wannier calculator

export ASE_PW2WANNIER_COMMAND="/path/to/pw2wannier PREFIX.pw2wan"

"""


import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'

class PW2Wannier(FileIOCalculator):
    """
    """
    implemented_properties = []
    command = 'pw2wannierPREFIX.pw2wan'

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
        io.write(self.label + '.pw2wan', atoms, **self.parameters)

    def read_results(self):
        return
        # output = io.read(self.label + '.pwo')
        # self.calc = output.calc
        # self.results = output.calc.results
