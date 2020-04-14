"""Module to read and write atoms in xtl file format for the muSTEM software.

See http://tcmp.ph.unimelb.edu.au/mustem/muSTEM.html for a few examples of
this format and the documentation of muSTEM.

See https://github.com/HamishGBrown/MuSTEM for the source code of muSTEM.
"""

import numpy as np

from ase.atoms import symbols2numbers
from ase.utils import reader


@reader
def read_mustem(fd):
    """Import muSTEM input file.

    Reads cell, atom positions, etc. from muSTEM xtl file
    """

    from ase import Atoms
    from ase.geometry import cellpar_to_cell

    # Read comment:
    fd.readline()

    # Parse unit cell parameter
    cellpar = [float(i) for i in fd.readline().split()[:3]]
    cell = cellpar_to_cell(cellpar)

    # beam energy
    fd.readline()

    # Number of different type of atoms
    element_number = int(fd.readline().strip())

    symbols = []
    positions = []

    for i in range(element_number):
        # Read the element
        symbol = str(fd.readline().strip())
        atoms_number = int(fd.readline().split()[0])
        # read all the position for each element
        for j in range(atoms_number):
            line = fd.readline()
            positions.append([float(i) for i in line.split()])
            symbols.append(symbol)

    atoms = Atoms(cell=cell, scaled_positions=positions)
    atoms.set_chemical_symbols(symbols)

    return atoms


class XtlmuSTEMWriter:
    """See the docstring of the `write_mustem` function.
    """

    def __init__(self, atoms, keV, DW=None, comment=None, occupancy=1.0,
                 fit_cell_to_atoms=False):
        cell = atoms.get_cell()
        if not cell.orthorhombic:
            raise ValueError('To export to this format, the cell needs to be '
                             'orthorhombic.')
        if cell.rank < 3:
            raise ValueError('To export to this format, the cell size needs '
                             'to be set: current cell is {}.'.format(cell))
        self.atoms = atoms.copy()
        self.atom_types = list(set(atoms.symbols))
        self.keV = keV
        self.comment = comment
        self.occupancy = self._get_occupancy(occupancy)
        self.DW = self._get_DW(DW)

        self.numbers = symbols2numbers(self.atom_types)
        if fit_cell_to_atoms:
            self.atoms.translate(-self.atoms.positions.min(axis=0))
            self.atoms.set_cell(self.atoms.positions.max(axis=0))

    def _get_occupancy(self, occupancy):
        if np.isscalar(occupancy):
            occupancy = {atom: occupancy for atom in self.atom_types}
        elif isinstance(occupancy, dict):
            self._check_key_dictionary(occupancy, 'occupancy')

        return occupancy

    def _get_DW(self, DW):
        if np.isscalar(DW):
            if len(self.atom_types) > 1:
                raise ValueError('This cell contains more then one type of '
                                 'atoms and the Debye-Waller factor needs to '
                                 'be provided for each atom using a '
                                 'dictionary.')
            DW = {self.atom_types[0]: DW}
        elif isinstance(DW, dict):
            self._check_key_dictionary(DW, 'DW')

        if DW is None:
            raise ValueError('Missing Debye-Waller factors. It can be '
                             'provided as a dictionary with symbols as key or '
                             'if the cell contains only a single type of '
                             'element, the Debye-Waller factor can also be '
                             'provided as float.')

        return DW

    def _check_key_dictionary(self, d, dict_name):
        # Check if we have enough key
        for key in self.atom_types:
            if key not in d:
                raise ValueError('Missing the {0} key in the `{1}` dictionary.'
                                 ''.format(key, dict_name))

    def _get_position_array_single_atom_type(self, number):
        # Get the scaled (reduced) position for a single atom type
        return self.atoms.get_scaled_positions()[np.where(
            self.atoms.numbers == number)]

    def _get_file_header(self):
        # 1st line: comment line
        if self.comment is None:
            s = "{0} atoms with chemical formula: {1}\n".format(
                len(self.atoms),
                self.atoms.get_chemical_formula())
        else:
            s = self.comment
        # 2nd line: lattice parameter
        s += "{} {} {} {} {} {}\n".format(
            *self.atoms.get_cell_lengths_and_angles().tolist())
        # 3td line: acceleration voltage
        s += "{}\n".format(self.keV)
        # 4th line: number of different atom
        s += "{}\n".format(len(self.atom_types))
        return s

    def _get_element_header(self, atom_type, number, atom_type_number,
                            occupancy, DW):
        return "{0}\n{1} {2} {3} {4}\n".format(atom_type, number,
                                               atom_type_number, occupancy, DW)

    def _get_file_end(self):
        return "Orientation\n   1 0 0\n   0 1 0\n   0 0 1\n"

    def write_to_file(self, f):
        if isinstance(f, str):
            f = open(f, 'w')

        f.write(self._get_file_header())
        for atom_type, number, occupancy in zip(self.atom_types,
                                                self.numbers,
                                                self.occupancy):
            positions = self._get_position_array_single_atom_type(number)
            atom_type_number = positions.shape[0]
            f.write(self._get_element_header(atom_type, atom_type_number,
                                             number, self.occupancy[atom_type],
                                             self.DW[atom_type]))
            for pos in positions:
                f.write('{0} {1} {2}\n'.format(pos[0], pos[1], pos[2]))

        f.write(self._get_file_end())


def write_mustem(filename, *args, **kwargs):
    """Write muSTEM input file.

    Parameters:

    atoms: Atoms object

    keV: float
        Energy of the electron beam in keV required for the image simulation.

    DW: float or dictionary of float with atom type as key
        Debye-Waller factor of each atoms.

    occupancy: float or dictionary of float with atom type as key (optional)
        Occupancy of each atoms. Default value is `1.0`.

    comment: str (optional)
        Comments to be written in the first line of the file. If not
        provided, write the total number of atoms and the chemical formula.

    fit_cell_to_atoms: bool (optional)
        If `True`, fit the cell to the atoms positions. If negative coordinates
        are present in the cell, the atoms are translated, so that all
        positions are positive. If `False` (default), the atoms positions and
        the cell are unchanged.
    """

    writer = XtlmuSTEMWriter(*args, **kwargs)
    writer.write_to_file(filename)
