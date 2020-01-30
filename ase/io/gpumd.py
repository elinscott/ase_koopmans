import numpy as np
from itertools import product
from ..neighborlist import NeighborList


def write_atoms_gpumd(fname, atoms, maximum_neighbors=None, cutoff=None,
                      velocities=None, groupings=None, use_triclinic=False):
    """
    Writes atoms into GPUMD input format.

    Parameters
    ----------
    fname : file | str
        Name of file to which to write atoms object
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
    symbols = atoms.get_chemical_symbols()
    types = {s: i for i, s in enumerate(set(symbols))}

    # atoms lines
    for a, atm in enumerate(atoms):
        t = types[atm.symbol]
        line = (' {}' * 5)[1:].format(t, *atm.position, atm.mass)
        if groupings is not None:
            for grouping in groupings:
                for i, group in enumerate(grouping):
                    if a in group:
                        line += ' {}'.format(i)
                        break
        lines.append(line)

    with open(fname, 'w') as f:
        all_lines = '\n'.join(lines)
        f.write(all_lines)

