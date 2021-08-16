"""Generic Quantum ESPRESSO Calculator

"""


import numpy as np
import warnings
from ase import io
from ase.spectrum.band_structure import BandStructure
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

    def get_number_of_spins(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Number of spins')
        nspins = self.calc.get_number_of_spins()
        if nspins is None:
            warnings.warn(warn_template % 'Number of spins')
        return nspins


class EspressoWithBandstructure:
    # Parent class to be used for Espresso calculators that yield bandstructures. These classes implement a
    # self.band_structure() which is called after self.calculate(). This is done because after self.calculate()
    # because we require access to the band path

    # Putting the band structure in self.results is very un-ASE-y, so we should ultimately align all of this
    # with the more general self.band_structure() method of ASE

    def eigenvalues_from_results(self):
        # Method that should return np array of eigenvalues
        raise NotImplementedError

    @property
    def vbm_index(self):
        raise NotImplementedError

    @property
    def vbm_energy(self):
        eigenvalues_np = self.eigenvalues_from_results()
        return np.max(eigenvalues_np[:, :, self.vbm_index])

    def band_structure(self, vbm_to_zero=True):
        # Fetch bandstructure from results
        eigenvalues_np = self.eigenvalues_from_results()

        if vbm_to_zero:
            # Shift so that the VBM = 0.0
            eigenvalues_np -= self.vbm_energy

        self.results['band structure'] = BandStructure(self.parameters['kpts'], eigenvalues_np)
        return self.results['band structure']
