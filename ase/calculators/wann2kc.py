""" Calculator for wann_to_kc.x, a code for performing Koopmans calculations in periodic systems

export ASE_KOOPMANS_HAM_COMMAND="/path/to/wann_to_kc.x -in PREFIX.khi > PREFIX.kho"

N.B. the extensions must be .w2ki and .w2ko

Run wann_to_kc.x jobs.
"""


import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'


class Wann2KC(FileIOCalculator):
    """
    """
    implemented_properties = []

    # Default command does not use parallelism and assumes wann_to_kc.x is on the user's path
    command = 'wann_to_kc.x -in PREFIX.w2ki > PREFIX.w2ko'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='wann2kc', atoms=None, **kwargs):
        """
        All options for wann_to_kc.x are copied verbatim to the input file, and put
        into the correct section. Use ``input_data`` for parameters that are
        already in a dict, all other ``kwargs`` are passed as parameters.

        input_data: dict
            A flat or nested dictionary with input parameters for wann_to_kc.x
        pseudopotentials: dict
            A filename for each atomic species, e.g.
            ``{'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}``.
            A dummy name will be used if none are given.
        kspacing: float
            Generate a grid of k-points with this as the minimum distance,
            in A^-1 between them in reciprocal space. If set to None, kpts
            will be used instead.
        kpts: (int, int, int), dict, or BandPath
            If kpts is a tuple (or list) of 3 integers, it is interpreted
            as the dimensions of a Monkhorst-Pack grid.
            If kpts is a dict, it will either be interpreted as a path
            in the Brillouin zone (*) if it contains the 'path' keyword,
            otherwise it is converted to a Monkhorst-Pack grid (**).
            (*) see ase.dft.kpoints.bandpath
            (**) see ase.calculators.calculator.kpts2sizeandoffsets
        koffset: (int, int, int)
            Offset of kpoints in each direction. Must be 0 (no offset) or
            1 (half grid offset). Setting to True is equivalent to (1, 1, 1).

        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if atoms is not None:
            self.atoms = atoms
            self.atoms.calc = self

        self.calc = None

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        io.write(self.label + '.w2ki', atoms, **self.parameters)

    def read_results(self):
        output = io.read(self.label + '.w2ko')
        self.calc = output.calc
        self.results = output.calc.results
