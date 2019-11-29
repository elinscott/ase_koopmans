"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os

from ase import io
from ase.calculators.calculator import FileIOCalculator


class NWChem(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'stress', 'dipole']
    command = 'nwchem PREFIX.nwi > PREFIX.nwo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='nwchem', atoms=None, command=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command, **kwargs)
        self.calc = None

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        # Prepare perm and scratch directories
        label = self.parameters.get('label', 'nwchem')
        perm = os.path.abspath(self.parameters.get('perm', label))
        scratch = os.path.abspath(self.parameters.get('scratch', label))
        os.makedirs(perm, exist_ok=True)
        os.makedirs(scratch, exist_ok=True)

        io.write(self.label + '.nwi', atoms, properties=properties,
                 label=self.label, **self.parameters)

    def read_results(self):
        output = io.read(self.label + '.nwo')
        self.calc = output.calc
        self.results = output.calc.results
