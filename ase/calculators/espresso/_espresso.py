"""Generic Quantum ESPRESSO Calculator

"""


import numpy as np
import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'


class EspressoParent(FileIOCalculator):
    """
    Parent class for all of the Quantum ESPRESSO calculators
    """
    implemented_properties = []
    ext_in = ''
    ext_out = ''

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='espresso', atoms=None, **kwargs):

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if atoms is not None:
            self.atoms = atoms
            self.atoms.calc = self

        self.calc = None

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.results = {}

    def write_input(self, atoms, properties=None, system_changes=None):
        # Create the appropriate directory
        super().write_input(atoms, properties, system_changes)
        # Write the input file
        io.write(self.label + self.ext_in, atoms, **self.parameters)

    def read_results(self):
        output = io.read(self.label + self.ext_out)
        self.calc = output.calc
        self.results = output.calc.results
        if hasattr(output.calc, 'kpts'):
            self.kpts = output.calc.kpts

    def get_number_of_spins(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Number of spins')
        nspins = self.calc.get_number_of_spins()
        if nspins is None:
            warnings.warn(warn_template % 'Number of spins')
        return nspins
