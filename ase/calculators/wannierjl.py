"""
WannierJL calculator
"""

import toml
from pathlib import Path

from ase import io
from ase.calculators.calculator import FileIOCalculator


class WannierJL(FileIOCalculator):
    implemented_properties = []
    command = 'wjl TASK --config=PREFIX.wjli PREFIX > PREFIX.wjlo 2>&1'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='wannier', atoms=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file, label, atoms, **kwargs)

        if atoms is not None:
            self.atoms = atoms
            self.atoms.calc = self

        self.calc = None

    def calculate(self, *args, **kwargs):
        # wjl "task" is specified via the command-line rather than an input file
        task = self.parameters.get('task', 'splitvc')
        self.command = self.command.replace('TASK', task)

        return super().calculate(*args, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):

        # Create the appropriate directory
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        # Assert that the W90 input file already exists
        w90_input_file = Path(self.label + '.win')
        if not w90_input_file.exists():
            raise FileNotFoundError(f'{w90_input_file} must exist before calling WannierJL')
    
        # Write the input file
        io.write(self.label + '.wjli', atoms)
        
    def read_results(self):
        output = io.read(self.label + '.wjlo')
        self.calc = output.calc
        self.results = output.calc.results
