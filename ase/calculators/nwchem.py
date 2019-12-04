"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os
import numpy as np

from ase import io
from ase.units import Hartree
from ase.calculators.calculator import FileIOCalculator
from ase.dft.band_structure import BandStructure


class NWChem(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'stress', 'dipole']
    command = 'nwchem PREFIX.nwi > PREFIX.nwo'
    accepts_bandpath_keyword = True

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
        perm = os.path.abspath(self.parameters.get('perm', self.label))
        scratch = os.path.abspath(self.parameters.get('scratch', self.label))
        os.makedirs(perm, exist_ok=True)
        os.makedirs(scratch, exist_ok=True)

        io.write(self.label + '.nwi', atoms, properties=properties,
                 label=self.label, **self.parameters)

    def read_results(self):
        output = io.read(self.label + '.nwo')
        self.calc = output.calc
        self.results = output.calc.results

    def band_structure(self):
        self.calculate()
        perm = self.parameters.get('perm', self.label)
        if self.calc.get_spin_polarized():
            alpha = np.loadtxt(os.path.join(perm, self.label + '.alpha_band'))
            beta = np.loadtxt(os.path.join(perm, self.label + '.beta_band'))
            energies = np.array([alpha[:, 1:], beta[:, 1:]]) * Hartree
        else:
            data = np.loadtxt(os.path.join(perm, self.label + '.restricted_band'))
            energies = data[np.newaxis, :, 1:] * Hartree
        eref = self.calc.get_fermi_level()
        if eref is None:
            eref = 0.
        return BandStructure(self.parameters.bandpath, energies, eref)
