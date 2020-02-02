import numpy as np
from itertools import product
from ..neighborlist import NeighborList
from ..data import atomic_masses, chemical_symbols


def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_nearest_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def write_gpumd(fileobj, atoms, maximum_neighbors=None, cutoff=None,
                      velocities=None, groupings=None, use_triclinic=False):
    """
    Writes atoms into GPUMD input format.

    Parameters
    ----------
    fileobj : file | str
        File object or name of file to which to write atoms object
    atoms : Atoms
        Input structure
    maximum_neighbors: int
        Maximum number of neighbors any atom can ever have (not relevant when
        using force constant potentials)
    cutoff: float
        Initial cutoff distance used for building the neighbor list (not
        relevant when using force constant potentials)
    velocities: list[list[float]]
        Initial velocity components for all atoms
    groupings : list[list[list[int]]]
        Groups into which the individual atoms should be divided in the form of
        a list of list of lists. Specifically, the outer list corresponds to
        the grouping methods, of which there can be three at the most, which
        contains a list of groups in the form of lists of site indices. The
        sum of the lengths of the latter must be the same as the total number
        of atoms.
    use_triclinic: bool
        Use format for triclinic cells

    Raises 
    ------
    ValueError
        Raised if parameters are incompatible
    """

    # Check velocties parameter
    if velocities is None:
        has_velocity = 0
    else:
        has_velocity = 1
        if len(velocities) != len(atoms)):
            raise ValueError('The number of velocities ({}) does not match the'
                             ' number of atoms'
                             ' ({})!'.format(len(velocities), len(atoms))

    # Check groupings parameter
    if groupings is None:
        number_of_grouping_methods = 0
    else:
        number_of_grouping_methods = len(groupings)
        if number_of_grouping_methods > 3
        raise ValueError('There can be no more than 3 grouping methods!')
        for g, grouping in enumerate(groupings):
            all_indices = [i for group in grouping for i in group]
            if len(all_indices) != len(atoms) or\
                    if set(all_indices) != set(range(len(atoms)):
                raise ValueError('The indices listed in grouping method {} are'
                                 ' not compatible with the input'
                                 ' structure!'.format(g)
    # try to estimate a good maximum_neighbors
    if maximum_neighbors is None:
        if cutoffs is None:
            cutoffs = 0.1
            maximum_neighbors = 1
        else:
            nl = NeighborList([cutoff/2]*len(atoms), skin=2, bothways=True)
            nl.update(atoms)
            maximum_neighbors = 0
            for atom in atoms:
                maximum_neighbors = max(maximum_neighbors,
                                        len(nl.get_neighbors(atom.index)[0]))
                maximum_neighbors *= 2

    # header with cell
    lines = []
    if atoms.cell.orthorhombic and not use_triclinic:
        head_lines = ['{} {} {} 0 {} {}'.format(len(atoms), maximum_neighbors,
                                                 cutoff, has_velocity,
                                                 number_of_grouping_methods)]
        head_lines.append((' {}' * 6)[1:].format(*atoms.pbc.astype(int),
                                                 *atoms.cell.lengths()))
    else:
        head_lines = ['{} {} {} 1 {} {}'.format(len(atoms), maximum_neighbors,
                                                 cutoff, has_velocity,
                                                 number_of_grouping_methods)]
        head_lines.append((' {}' * 12)[1:].format(*atoms.pbc.astype(int),
                                                  *atoms.cell[:].flatten()))
    lines += head_lines

    # symbols to integers starting at 0
    symbol_type_map = {}
    for symbol in atoms.get_chemical_symbols():
        if symbol not in symbol_type_map:
            symbol_type_map[symbol] = len(symbol_type_map)

    # atoms lines
    for a, atm in enumerate(atoms):
        t = symbol_type_map[atm.symbol]
        line = (' {}' * 5)[1:].format(t, *atm.position, atm.mass)
        if groupings is not None:
            for grouping in groupings:
                for i, group in enumerate(grouping):
                    if a in group:
                        line += ' {}'.format(i)
                        break
        lines.append(line)

    with open(fileobj, 'w') as f:
        all_lines = '\n'.join(lines)
        f.write(all_lines)


def read_xyz_input_gpumd(fileobj='xyz.in', species_types=None,
                         isotope_masses=None):

    """
    Read the structure input file for GPUMD and return an ase Atoms object
    togehter with a dictionary with parameters and a types-to-symbols map

    Parameters
    ----------
    fileobj : file | str
        File object or name of file from which to read the Atoms object
    species_types : List[str]
        List with the chemical symbols that correspond to each type, will take
        precedence over isotope_masses
    isotope_masses: Dict[str, List[float]]
        Dictionary with chemical symbols and lists of the associated atomic
        masses, which is used to identify the chemical symbols that correspond
        to the types not found in species_types. The default is to find the 
        closest match :data:`ase.data.atomic_masses`.

    Returns
    -------
    atoms : Atoms
        Atoms object
    input_parameters : Dict[str, int]
        Dictionary with parameters from the first row of the input file, namely
        'N', 'M', 'cutoff', 'use_triclinic', 'has_velocity' and 'num_of_groups'
    type_symbol_map : Dict[int, str]
        Dictionary with types and the corresponding chemical symbols

    Raises 
    ------
    ValueError
        Raised if the list of species is incompatible with the input file
    """
    if isotope_masses is not None:
        mass_symbols = {mass: symbol for symbol, masses in
                        isotope_masses.items() for mass in masses}

    # Read file
    with open(fileobj) as f:
        first_line = np.loadtxt(f, max_rows=1)
        second_line = np.loadtxt(f, max_rows=1)
        xyz = np.loadtxt(f)

    # Parse first line
    input_parameters = {}
    keys = ['N', 'M', 'cutoff', 'use_triclinic', 'has_velocity',
            'num_of_groups']
    types = [float if key == 'cutoff' else int for key in keys]
    for k, (key, typ) in enumerate(zip(keys, types)):
        input_parameters[key] = typ(first_line[k])

    # Parse second line
    pbc = second_line[:3].astype(int)
    if input_parameters['use_triclinic']:
        cell = second_line[3:].reshape((3, 3))
    else:
        cell = np.diag(second_line[3:])

    # Initiate the ase Atoms object
    info = dict()
    atoms = Atoms()
    atoms.set_pbc(pbc)
    atoms.set_cell(cell)

    if species_types is not None:
        if len(species_types) > input_parameters['N']:
             raise ValueError('The number of species types ({}) exceeds the'
                              ' number of atom types'
                              ' {}'.format(len(species_types),
                                           input_parameters['N'])
        type_symbol_map = {index: symbol for index, symbol in
                            enumerate(species_types)}
    else:
        type_symbol_map = {}
    for i, xyz_row in enumerate(xyz):
        # Determine the atomic species from the mass
        atom_type = xyz_row[0]
        mass = xyz_row[4]
        if atom_type not in type_symbol_map:
            if isotope_masses is not None:
                nearest_mass = find_nearest_value(mass_symbols.keys(), mass)
                symbol = mass_symbols[nearest_mass]
            else:
                symbol = chemical_symbols[
                    find_nearest_index(atomic_masses, mass)]
            type_symbol_map[atom_type] = symbol

        # Create and add an ase Atom object
        position = xyz_row[1:4]
        symbol = type_symbol_map[atom_type]
        atom = Atom(symbol, position, mass=mass)
        atoms.append(atom)

        # Collect data regarding velocities and groups
        data = dict()
        if input_parameters['has_velocity']:
            data['velocity'] = xyz_row[5:8]
        if input_parameters['num_of_groups']:
            data['groups'] = xyz_row[8:].astype(int)
        info[i] = data

    # Add data regarding velocities and groups
    atoms.info = info

    return atoms, input_parameters, type_symbol_map


def read_gpumd(fileobj='xyz.in', species_types=None, isotope_masses=None):
    """
    Read Atoms object from a GPUMD structure input file

    Parameters
    ----------
    fileobj : file | str
        File object or name of file from which to read the Atoms object
    species_types : List[str]
        List with the chemical symbols that correspond to each type, will take
        precedence over isotope_masses
    isotope_masses: Dict[str, List[float]]
        Dictionary with chemical symbols and lists of the associated atomic
        masses, which is used to identify the chemical symbols that correspond
        to the types not found in species_types. The default is to find the 
        closest match :data:`ase.data.atomic_masses`.

    Returns
    -------
    atoms : Atoms
        Atoms object

    Raises 
    ------
    ValueError
        Raised if the list of species is incompatible with the input file
    """
   
    return read_xyz_input_gpumd(fileobj, species_types, isotope_masses)[0]
