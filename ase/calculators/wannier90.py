""" wannier90 calculator

export ASE_WANNIER90_COMMAND="/path/to/wannier90 PREFIX.win"

"""

import os
import numpy as np
from ase import io
from ase.dft.kpoints import bandpath, BandPath
from ase.spectrum.band_structure import BandStructure
from ase.calculators.calculator import FileIOCalculator


class Wannier90(FileIOCalculator):
    """
    """
    implemented_properties = []
    command = 'wannier90.x PREFIX.win > PREFIX.wout 2>&1'

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

        if self.parameters.get('bands_plot', False):
            self.read_bandstructure()

    def read_bandstructure(self):
        if not os.path.isfile(f'{self.label}_band.dat'):
            return

        # Construct the higher-resolution bandpath from the *_band.labelinfo.dat file
        with open(f'{self.label}_band.labelinfo.dat') as fd:
            flines = fd.readlines()
        kpts = []

        self.atoms.cell.pbc = True
        for start, end in zip(flines[:-1], flines[1:]):
            start_label, i_start = start.split()[:2]
            end_label, i_end = end.split()[:2]
            kpts += self.atoms.cell.bandpath(start_label + end_label, int(i_end) - int(i_start) + 1).kpts[:-1].tolist()
        kpts.append([float(x) for x in flines[-1].split()[-3:]])
        path = self.parameters.kpoint_path.path
        special_points = self.parameters.kpoint_path.special_points
        kpath = BandPath(self.atoms.cell, kpts, path=path, special_points=special_points)

        # Read in the eigenvalues from the *_band.dat file
        with open(f'{self.label}_band.dat') as fd:
            flines = fd.readlines()
        eigs = [[]]
        for line in flines[:-1]:
            splitline = line.strip().split()
            if len(splitline) == 0:
                eigs.append([])
            else:
                eigs[-1].append(line.strip().split()[-1])
        eigs = np.array(eigs, dtype=float).T

        # Construct the bandstructure
        bs = BandStructure(kpath, eigs[np.newaxis, :, :])

        self.results['band structure'] = bs
