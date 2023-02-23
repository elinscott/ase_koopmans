""" wannier90 calculator

export ASE_WANNIERJL_COMMAND="/path/to/wjl TASK PREFIX"

"""

import os
import numpy as np
from ase import io
from ase.dft.kpoints import bandpath, BandPath
from ase.spectrum.band_structure import BandStructure
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.wannier90 import write_input as write_w90_input


class WannierJL(FileIOCalculator):
    """
    """
    implemented_properties = []
    command = 'wjl TASK PREFIX > PREFIX.wjlo 2>&1'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='wannier', atoms=None, **kwargs):
        """
        All options for wannier90 are copied verbatim to the input file

        """
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if atoms is not None:
            self.atoms = atoms
            self.atoms.calc = self

        self.calc = None

    def calculate(self, *args, **kwargs):
        task = self.parameters.pop('task', 'splitvc')
        self.command = self.command.replace('TASK', task)
        return super().calculate(*args, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        # Create the appropriate directory
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        # Write the W90 input file
        io.write(self.label + '.win', atoms)

    def read_results(self):
        output = io.read(self.label + '.wjlo')
        self.calc = output.calc
        self.results = output.calc.results
