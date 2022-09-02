""" Utilities for Quantum ESPRESSO file I/O

Written by Edward Linscott 2020-21

"""

import warnings
import os
from ase.dft.kpoints import BandPath, labels_from_kpts
import numpy as np
from pathlib import Path
import operator as op
from collections import OrderedDict
from ase.units import create_units
from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
from ase.constraints import FixAtoms, FixCartesian
from ase.data import atomic_numbers, chemical_symbols

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')


class Namelist(OrderedDict):
    """Case insensitive dict that emulates Fortran Namelists."""

    def __contains__(self, key):
        return super(Namelist, self).__contains__(key.lower())

    def __delitem__(self, key):
        return super(Namelist, self).__delitem__(key.lower())

    def __getitem__(self, key):
        return super(Namelist, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        super(Namelist, self).__setitem__(key.lower(), value)

    def get(self, key, default=None):
        return super(Namelist, self).get(key.lower(), default)


# Number of valence electrons in the pseudopotentials recommended by
# http://materialscloud.org/sssp/. These are just used as a fallback for
# calculating initial magetization values which are given as a fraction
# of valence electrons.
SSSP_VALENCE = [
    0, 1.0, 2.0, 3.0, 4.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 3.0, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
    18.0, 19.0, 20.0, 13.0, 14.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 12.0, 13.0, 14.0, 15.0, 6.0,
    7.0, 18.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
    19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 36.0, 27.0, 14.0, 15.0, 30.0,
    15.0, 32.0, 19.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0]


class Namelist(OrderedDict):
    """Case insensitive dict that emulates Fortran Namelists."""

    def __contains__(self, key):
        return super(Namelist, self).__contains__(key.lower())

    def __delitem__(self, key):
        return super(Namelist, self).__delitem__(key.lower())

    def __getitem__(self, key):
        return super(Namelist, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        super(Namelist, self).__setitem__(key.lower(), value)

    def get(self, key, default=None):
        return super(Namelist, self).get(key.lower(), default)


def get_kpoints(card_lines, cell=None):
    # Find kpts, koffset from card_lines
    for i, line in enumerate(card_lines):
        if line.startswith('K_POINTS'):
            mode = line.split()[-1].strip('{}')
            if mode.lower() == 'gamma':
                return [1, 1, 1], [0, 0, 0], True
            elif mode.lower() == 'automatic':
                splitline = card_lines[i + 1].strip().split()
                kpts = [int(x) for x in splitline[:3]]
                koffset = [int(x) for x in splitline[3:]]
                return kpts, koffset, False
            elif mode.lower() == 'crystal_b':
                assert cell is not None
                n_kpts = int(card_lines[i + 1])
                kpts = np.array([[float(x) for x in line.split()[:3]] for line in card_lines[i + 2: i + 2 + n_kpts]])
                _, _, path_list = labels_from_kpts(kpts, cell)
                return BandPath(path=''.join(path_list), cell=cell, special_points=cell.bandpath().special_points,
                                kpts=kpts), [0, 0, 0], False
            else:
                raise ValueError('Failed to parse K_POINTS block')
    raise ValueError('Failed to find K_POINTS block')


def ibrav_to_cell(system):
    """
    Convert a value of ibrav to a cell. Any unspecified lattice dimension
    is set to 0.0, but will not necessarily raise an error. Also return the
    lattice parameter.

    Parameters
    ----------
    system : dict
        The &SYSTEM section of the input file, containing the 'ibrav' setting,
        and either celldm(1)..(6) or a, b, c, cosAB, cosAC, cosBC.

    Returns
    -------
    alat, cell : float, np.array
        Cell parameter in Angstrom, and
        The 3x3 array representation of the cell.

    Raises
    ------
    KeyError
        Raise an error if any required keys are missing.
    NotImplementedError
        Only a limited number of ibrav settings can be parsed. An error
        is raised if the ibrav interpretation is not implemented.
    """
    if 'celldm(1)' in system and 'a' in system:
        raise KeyError('do not specify both celldm and a,b,c!')
    elif 'celldm(1)' in system:
        # celldm(x) in bohr
        alat = system['celldm(1)'] * units['Bohr']
        b_over_a = system.get('celldm(2)', 0.0)
        c_over_a = system.get('celldm(3)', 0.0)
        cosab = system.get('celldm(4)', 0.0)
        cosac = system.get('celldm(5)', 0.0)
        cosbc = 0.0
        if system['ibrav'] == 14:
            cosbc = system.get('celldm(4)', 0.0)
            cosac = system.get('celldm(5)', 0.0)
            cosab = system.get('celldm(6)', 0.0)
    elif 'a' in system:
        # a, b, c, cosAB, cosAC, cosBC in Angstrom
        alat = system['a']
        b_over_a = system.get('b', 0.0) / alat
        c_over_a = system.get('c', 0.0) / alat
        cosab = system.get('cosab', 0.0)
        cosac = system.get('cosac', 0.0)
        cosbc = system.get('cosbc', 0.0)
    else:
        raise KeyError("Missing celldm(1) or a cell parameter.")

    if system['ibrav'] == 1:
        cell = np.identity(3) * alat
    elif system['ibrav'] == 2:
        cell = np.array([[-1.0, 0.0, 1.0],
                         [0.0, 1.0, 1.0],
                         [-1.0, 1.0, 0.0]]) * (alat / 2)
    elif system['ibrav'] == 3:
        cell = np.array([[1.0, 1.0, 1.0],
                         [-1.0, 1.0, 1.0],
                         [-1.0, -1.0, 1.0]]) * (alat / 2)
    elif system['ibrav'] == -3:
        cell = np.array([[-1.0, 1.0, 1.0],
                         [1.0, -1.0, 1.0],
                         [1.0, 1.0, -1.0]]) * (alat / 2)
    elif system['ibrav'] == 4:
        cell = np.array([[1.0, 0.0, 0.0],
                         [-0.5, 0.5 * 3**0.5, 0.0],
                         [0.0, 0.0, c_over_a]]) * alat
    elif system['ibrav'] == 5:
        tx = ((1.0 - cosab) / 2.0)**0.5
        ty = ((1.0 - cosab) / 6.0)**0.5
        tz = ((1 + 2 * cosab) / 3.0)**0.5
        cell = np.array([[tx, -ty, tz],
                         [0, 2 * ty, tz],
                         [-tx, -ty, tz]]) * alat
    elif system['ibrav'] == -5:
        ty = ((1.0 - cosab) / 6.0)**0.5
        tz = ((1 + 2 * cosab) / 3.0)**0.5
        a_prime = alat / 3**0.5
        u = tz - 2 * 2**0.5 * ty
        v = tz + 2**0.5 * ty
        cell = np.array([[u, v, v],
                         [v, u, v],
                         [v, v, u]]) * a_prime
    elif system['ibrav'] == 6:
        cell = np.array([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, c_over_a]]) * alat
    elif system['ibrav'] == 7:
        cell = np.array([[1.0, -1.0, c_over_a],
                         [1.0, 1.0, c_over_a],
                         [-1.0, -1.0, c_over_a]]) * (alat / 2)
    elif system['ibrav'] == 8:
        cell = np.array([[1.0, 0.0, 0.0],
                         [0.0, b_over_a, 0.0],
                         [0.0, 0.0, c_over_a]]) * alat
    elif system['ibrav'] == 9:
        cell = np.array([[1.0 / 2.0, b_over_a / 2.0, 0.0],
                         [-1.0 / 2.0, b_over_a / 2.0, 0.0],
                         [0.0, 0.0, c_over_a]]) * alat
    elif system['ibrav'] == -9:
        cell = np.array([[1.0 / 2.0, -b_over_a / 2.0, 0.0],
                         [1.0 / 2.0, b_over_a / 2.0, 0.0],
                         [0.0, 0.0, c_over_a]]) * alat
    elif system['ibrav'] == 10:
        cell = np.array([[1.0 / 2.0, 0.0, c_over_a / 2.0],
                         [1.0 / 2.0, b_over_a / 2.0, 0.0],
                         [0.0, b_over_a / 2.0, c_over_a / 2.0]]) * alat
    elif system['ibrav'] == 11:
        cell = np.array([[1.0 / 2.0, b_over_a / 2.0, c_over_a / 2.0],
                         [-1.0 / 2.0, b_over_a / 2.0, c_over_a / 2.0],
                         [-1.0 / 2.0, -b_over_a / 2.0, c_over_a / 2.0]]) * alat
    elif system['ibrav'] == 12:
        sinab = (1.0 - cosab**2)**0.5
        cell = np.array([[1.0, 0.0, 0.0],
                         [b_over_a * cosab, b_over_a * sinab, 0.0],
                         [0.0, 0.0, c_over_a]]) * alat
    elif system['ibrav'] == -12:
        sinac = (1.0 - cosac**2)**0.5
        cell = np.array([[1.0, 0.0, 0.0],
                         [0.0, b_over_a, 0.0],
                         [c_over_a * cosac, 0.0, c_over_a * sinac]]) * alat
    elif system['ibrav'] == 13:
        sinab = (1.0 - cosab**2)**0.5
        cell = np.array([[1.0 / 2.0, 0.0, -c_over_a / 2.0],
                         [b_over_a * cosab, b_over_a * sinab, 0.0],
                         [1.0 / 2.0, 0.0, c_over_a / 2.0]]) * alat
    elif system['ibrav'] == 14:
        sinab = (1.0 - cosab**2)**0.5
        v3 = [c_over_a * cosac,
              c_over_a * (cosbc - cosac * cosab) / sinab,
              c_over_a * ((1 + 2 * cosbc * cosac * cosab
                           - cosbc**2 - cosac**2 - cosab**2)**0.5) / sinab]
        cell = np.array([[1.0, 0.0, 0.0],
                         [b_over_a * cosab, b_over_a * sinab, 0.0],
                         v3]) * alat
    else:
        raise NotImplementedError('ibrav = {0} is not implemented'
                                  ''.format(system['ibrav']))

    return alat, cell


def get_atomic_positions(lines, n_atoms, cell=None, alat=None):
    """Parse atom positions from ATOMIC_POSITIONS card.

    Parameters
    ----------
    lines : list[str]
        A list of lines containing the ATOMIC_POSITIONS card.
    n_atoms : int
        Expected number of atoms. Only this many lines will be parsed.
    cell : np.array
        Unit cell of the crystal. Only used with crystal coordinates.
    alat : float
        Lattice parameter for atomic coordinates. Only used for alat case.

    Returns
    -------
    positions : list[(str, (float, float, float), (float, float, float))]
        A list of the ordered atomic positions in the format:
        label, (x, y, z), (if_x, if_y, if_z)
        Force multipliers are set to None if not present.

    Raises
    ------
    ValueError
        Any problems parsing the data result in ValueError

    """

    positions = None
    # no blanks or comment lines, can the consume n_atoms lines for positions
    trimmed_lines = (line for line in lines
                     if line.strip() and not line[0] == '#')

    for line in trimmed_lines:
        if line.strip().startswith('ATOMIC_POSITIONS'):
            if positions is not None:
                raise ValueError('Multiple ATOMIC_POSITIONS specified')
            # Priority and behaviour tested with QE 5.3
            if 'crystal_sg' in line.lower():
                raise NotImplementedError('CRYSTAL_SG not implemented')
            elif 'crystal' in line.lower():
                cell = cell
            elif 'bohr' in line.lower():
                cell = np.identity(3) * units['Bohr']
            elif 'angstrom' in line.lower():
                cell = np.identity(3)
            # elif 'alat' in line.lower():
            #     cell = np.identity(3) * alat
            else:
                if alat is None:
                    raise ValueError('Set lattice parameter in &SYSTEM for '
                                     'alat coordinates')
                # Always the default, will be DEPRECATED as mandatory
                # in future
                cell = np.identity(3) * alat

            positions = []
            for _dummy in range(n_atoms):
                split_line = next(trimmed_lines).split()
                # These can be fractions and other expressions
                position = np.dot((infix_float(split_line[1]),
                                   infix_float(split_line[2]),
                                   infix_float(split_line[3])), cell)
                if len(split_line) > 4:
                    force_mult = (float(split_line[4]),
                                  float(split_line[5]),
                                  float(split_line[6]))
                else:
                    force_mult = None

                positions.append((split_line[0], position, force_mult))

    return positions


def get_cell_parameters(lines, alat=None):
    """Parse unit cell from CELL_PARAMETERS card.

    Parameters
    ----------
    lines : list[str]
        A list with lines containing the CELL_PARAMETERS card.
    alat : float | None
        Unit of lattice vectors in Angstrom. Only used if the card is
        given in units of alat. alat must be None if CELL_PARAMETERS card
        is in Bohr or Angstrom. For output files, alat will be parsed from
        the card header and used in preference to this value.

    Returns
    -------
    cell : np.array | None
        Cell parameters as a 3x3 array in Angstrom. If no cell is found
        None will be returned instead.
    cell_alat : float | None
        If a value for alat is given in the card header, this is also
        returned, otherwise this will be None.

    Raises
    ------
    ValueError
        If CELL_PARAMETERS are given in units of bohr or angstrom
        and alat is not
    """

    cell = None
    cell_alat = None
    # no blanks or comment lines, can take three lines for cell
    trimmed_lines = (line for line in lines
                     if line.strip() and not line[0] == '#')

    for line in trimmed_lines:
        if line.strip().startswith('CELL_PARAMETERS'):
            if cell is not None:
                # multiple definitions
                raise ValueError('CELL_PARAMETERS specified multiple times')
            # Priority and behaviour tested with QE 5.3
            if 'bohr' in line.lower():
                if alat is not None:
                    raise ValueError('Lattice parameters given in '
                                     '&SYSTEM celldm/A and CELL_PARAMETERS '
                                     'bohr')
                cell_units = units['Bohr']
            elif 'angstrom' in line.lower():
                cell_units = 1.0
            elif 'alat' in line.lower():
                # Output file has (alat = value) (in Bohrs)
                if '=' in line:
                    alat = float(line.strip(') \n').split()[-1]) * units['Bohr']
                    cell_alat = alat
                elif alat is None:
                    raise ValueError('Lattice parameters must be set in '
                                     '&SYSTEM for alat units')
                cell_units = alat
            elif alat is None:
                # may be DEPRECATED in future
                cell_units = units['Bohr']
            else:
                # may be DEPRECATED in future
                cell_units = alat
            # Grab the parameters; blank lines have been removed
            cell = [[ffloat(x) for x in next(trimmed_lines).split()[:3]],
                    [ffloat(x) for x in next(trimmed_lines).split()[:3]],
                    [ffloat(x) for x in next(trimmed_lines).split()[:3]]]
            cell = np.array(cell) * cell_units

    return cell, cell_alat


def get_pseudopotentials(lines, n_types):
    """Parse atom positions from ATOMIC_SPECIES card.

    Parameters
    ----------
    lines : list[str]
        A list of lines containing the ATOMIC_SPECIES card.
    n_types : int
        Expected number of species. Only this many lines will be parsed.

    Returns
    -------
    pseudos : {str : str, ...}
        A dictionary of the species and corresponding pseudopotentials

    Raises
    ------
    ValueError
        Any problems parsing the data result in ValueError

    """

    i_start = [l.lower() for l in lines].index('atomic_species') + 1
    pseudos = {l.split()[0]: l.split()[2] for l in lines[i_start:i_start + n_types]}

    return pseudos


def str_to_value(string):
    """Attempt to convert string into int, float (including fortran double),
    or bool, in that order, otherwise return the string.
    Valid (case-insensitive) bool values are: '.true.', '.t.', 'true'
    and 't' (or false equivalents).

    Parameters
    ----------
    string : str
        Test to parse for a datatype

    Returns
    -------
    value : any
        Parsed string as the most appropriate datatype of int, float,
        bool or string.

    """

    # Just an integer
    try:
        return int(string)
    except ValueError:
        pass
    # Standard float
    try:
        return float(string)
    except ValueError:
        pass
    # Fortran double
    try:
        return ffloat(string)
    except ValueError:
        pass

    # possible bool, else just the raw string
    if string.lower() in ('.true.', '.t.', 'true', 't'):
        return True
    elif string.lower() in ('.false.', '.f.', 'false', 'f'):
        return False
    else:
        return string.strip("'")


def read_fortran_namelist(fileobj):
    """Takes a fortran-namelist formatted file and returns nested
    dictionaries of sections and key-value data, followed by a list
    of lines of text that do not fit the specifications.

    Behaviour is taken from Quantum ESPRESSO 5.3. Parses fairly
    convoluted files the same way that QE should, but may not get
    all the MANDATORY rules and edge cases for very non-standard files:
        Ignores anything after '!' in a namelist, split pairs on ','
        to include multiple key=values on a line, read values on section
        start and end lines, section terminating character, '/', can appear
        anywhere on a line.
        All of these are ignored if the value is in 'quotes'.

    Parameters
    ----------
    fileobj : file
        An open file-like object.

    Returns
    -------
    data : dict of dict
        Dictionary for each section in the namelist with key = value
        pairs of data.
    card_lines : list of str
        Any lines not used to create the data, assumed to belong to 'cards'
        in the input file.

    """
    # Espresso requires the correct order
    data = Namelist()
    card_lines = []
    in_namelist = False
    section = 'none'  # can't be in a section without changing this

    for line in fileobj:
        # leading and trailing whitespace never needed
        line = line.strip()
        if line.startswith('&'):
            # inside a namelist
            section = line.split()[0][1:].lower()  # case insensitive
            if section in data:
                # Repeated sections are completely ignored.
                # (Note that repeated keys overwrite within a section)
                section = "_ignored"
            data[section] = Namelist()
            in_namelist = True
        if not in_namelist and line:
            # Stripped line is Truthy, so safe to index first character
            if line[0] not in ('!', '#'):
                card_lines.append(line)
        if in_namelist:
            # parse k, v from line:
            key = []
            value = None
            in_quotes = False
            for character in line:
                if character == ',' and value is not None and not in_quotes:
                    # finished value:
                    data[section][''.join(key).strip()] = str_to_value(
                        ''.join(value).strip())
                    key = []
                    value = None
                elif character == '=' and value is None and not in_quotes:
                    # start writing value
                    value = []
                elif character == "'":
                    # only found in value anyway
                    in_quotes = not in_quotes
                    value.append("'")
                elif character == '!' and not in_quotes:
                    break
                elif character == '/' and not in_quotes:
                    in_namelist = False
                    break
                elif value is not None:
                    value.append(character)
                else:
                    key.append(character)
            if value is not None:
                data[section][''.join(key).strip()] = str_to_value(
                    ''.join(value).strip())

    return data, card_lines


def ffloat(string):
    """Parse float from fortran compatible float definitions.

    In fortran exponents can be defined with 'd' or 'q' to symbolise
    double or quad precision numbers. Double precision numbers are
    converted to python floats and quad precision values are interpreted
    as numpy longdouble values (platform specific precision).

    Parameters
    ----------
    string : str
        A string containing a number in fortran real format

    Returns
    -------
    value : float | np.longdouble
        Parsed value of the string.

    Raises
    ------
    ValueError
        Unable to parse a float value.

    """

    if 'q' in string.lower():
        return np.longdouble(string.lower().replace('q', 'e'))
    else:
        return float(string.lower().replace('d', 'e'))


def label_to_symbol_and_tag(label):
    """Convert a valid espresso ATOMIC_SPECIES label to a
    chemical symbol.

    Parameters
    ----------
    label : str
        chemical symbol X (1 or 2 characters, case-insensitive)
        or chemical symbol plus a number or a letter, as in
        "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h;
        max total length cannot exceed 3 characters).

    Returns
    -------
    symbol : str
        The best matching species from ase.utils.chemical_symbols

    tag : int
        The tag provided

    Raises
    ------
    KeyError
        Couldn't find an appropriate species.

    Notes
    -----
        It's impossible to tell whether e.g. He is helium
        or hydrogen labelled 'e'.

        This function enforces tagging with integers
    """
    # possibly a two character species
    # ase Atoms need proper case of chemical symbols.

    symbol = None
    tag = None
    if len(label) >= 2:
        test_symbol = label[0].upper() + label[1].lower()
        if test_symbol in chemical_symbols:
            symbol, tag = test_symbol, label[2:]

    if symbol is None:
        # finally try with one character
        test_symbol = label[0].upper()
        if test_symbol in chemical_symbols:
            symbol, tag = test_symbol, label[1:]

    if symbol is None:
        raise KeyError('Could not parse species from label {0}.'
                       ''.format(label))

    # Enforce tag to be a positive integer (we will internally use 0 for non-tagged atoms)
    if tag == '':
        tag = 0
    elif tag == '0':
        raise ValueError(f'Please do not use 0 as a tag for {symbol}; instead use positive integers')
    else:
        try:
            tag = int(tag)
        except ValueError:
            raise ValueError(f'Please use positive integers as tags, rather than {tag} for {symbol}')

    return symbol, tag


def label_to_symbol(label):
    symbol, _ = label_to_symbol_and_tag(label)
    return symbol


def label_to_tag(label):
    _, tag = label_to_symbol_and_tag(label)
    return tag


def infix_float(text):
    """Parse simple infix maths into a float for compatibility with
    Quantum ESPRESSO ATOMIC_POSITIONS cards. Note: this works with the
    example, and most simple expressions, but the capabilities of
    the two parsers are not identical. Will also parse a normal float
    value properly, but slowly.

    >>> infix_float('1/2*3^(-1/2)')
    0.28867513459481287

    Parameters
    ----------
    text : str
        An arithmetic expression using +, -, *, / and ^, including brackets.

    Returns
    -------
    value : float
        Result of the mathematical expression.

    """

    def middle_brackets(full_text):
        """Extract text from innermost brackets."""
        start, end = 0, len(full_text)
        for (idx, char) in enumerate(full_text):
            if char == '(':
                start = idx
            if char == ')':
                end = idx + 1
                break
        return full_text[start:end]

    def eval_no_bracket_expr(full_text):
        """Calculate value of a mathematical expression, no brackets."""
        exprs = [('+', op.add), ('*', op.mul),
                 ('/', op.truediv), ('^', op.pow)]
        full_text = full_text.lstrip('(').rstrip(')')
        try:
            return float(full_text)
        except ValueError:
            for symbol, func in exprs:
                if symbol in full_text:
                    left, right = full_text.split(symbol, 1)  # single split
                    return func(eval_no_bracket_expr(left),
                                eval_no_bracket_expr(right))

    while '(' in text:
        middle = middle_brackets(text)
        text = text.replace(middle, '{}'.format(eval_no_bracket_expr(middle)))

    return float(eval_no_bracket_expr(text))

###
# Input file writing
###


def generic_construct_namelist(parameters=None, warn=False, local_keys=None, **kwargs):
    """
    Construct an ordered Namelist containing all the parameters given (as
    a dictionary or kwargs). Keys will be inserted into their appropriate
    section in the namelist and the dictionary may contain flat and nested
    structures. Any kwargs that match input keys will be incorporated into
    their correct section. All matches are case-insensitive, and returned
    Namelist object is a case-insensitive dict.

    If a key is not known to ase, but in a section within `parameters`,
    it will be assumed that it was put there on purpose and included
    in the output namelist. Anything not in a section will be ignored (set
    `warn` to True to see ignored keys).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    The priority of the keys is:
        kwargs[key] > parameters[key] > parameters[section][key]
    Only the highest priority item will be included.

    Parameters
    ----------
    parameters: dict
        Flat or nested set of input parameters.
    warn: bool
        Enable warnings for unused keys.

    Returns
    -------
    input_namelist: Namelist
        pw.x compatible namelist of input parameters.

    """
    # Convert everything to Namelist early to make case-insensitive
    if parameters is None:
        parameters = Namelist()
    else:
        # Maximum one level of nested dict
        # Don't modify in place
        parameters_namelist = Namelist()
        for key, value in parameters.items():
            if isinstance(value, dict):
                parameters_namelist[key] = Namelist(value)
            else:
                parameters_namelist[key] = value
        parameters = parameters_namelist

    # Just a dict
    kwargs = Namelist(kwargs)

    # Final parameter set
    input_namelist = Namelist()

    # Collect
    for section in local_keys:
        sec_list = Namelist()
        for key in local_keys[section]:
            # Check all three separately and pop them all so that
            # we can check for missing values later

            value = None
            if key in parameters.get(section, {}):
                value = parameters[section].pop(key)
            if key in parameters:
                value = parameters.pop(key)
            if key in kwargs:
                value = kwargs.pop(key)

            if value is not None:
                if isinstance(value, Path):
                    value = str(value) + os.path.sep
                sec_list[key] = value

            # Check if there is a key(i) version (no extra parsing)
            cp_parameters = parameters.copy()
            for arg_key in cp_parameters:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = parameters.pop(arg_key)
            cp_kwargs = kwargs.copy()
            for arg_key in cp_kwargs:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = kwargs.pop(arg_key)

        # Add to output
        input_namelist[section] = sec_list

    unused_keys = list(kwargs)
    # pass anything else already in a section
    for key, value in parameters.items():
        if key in local_keys and isinstance(value, dict):
            input_namelist[key].update(value)
        elif isinstance(value, dict):
            unused_keys.extend(list(value))
        else:
            unused_keys.append(key)

    if warn and unused_keys:
        warnings.warn('Unused keys: {}'.format(', '.join(unused_keys)))

    return input_namelist


def grep_valence(pseudopotential):
    """
    Given a UPF pseudopotential file, find the number of valence atoms.

    Parameters
    ----------
    pseudopotential: str
        Filename of the pseudopotential.

    Returns
    -------
    valence: float
        Valence as reported in the pseudopotential.

    Raises
    ------
    ValueError
        If valence cannot be found in the pseudopotential.
    """

    # Example lines
    # Sr.pbe-spn-rrkjus_psl.1.0.0.UPF:        z_valence="1.000000000000000E+001"
    # C.pbe-n-kjpaw_psl.1.0.0.UPF (new ld1.x):
    #                            ...PBC" z_valence="4.000000000000e0" total_p...
    # C_ONCV_PBE-1.0.upf:                     z_valence="    4.00"
    # Ta_pbe_v1.uspp.F.UPF:   13.00000000000      Z valence

    with open(pseudopotential) as psfile:
        for line in psfile:
            if 'z valence' in line.lower():
                return float(line.split()[0])
            elif 'z_valence' in line.lower():
                if line.split()[0] == '<PP_HEADER':
                    line = list(filter(lambda x: 'z_valence' in x,
                                       line.split(' ')))[0]
                return float(line.split('=')[-1].strip().strip('"'))
        else:
            raise ValueError('Valence missing in {}'.format(pseudopotential))


def cell_to_ibrav(cell, ibrav):
    """
    Calculate the appropriate `celldm(..)` parameters for the given ibrav
    using the given cell. The units for `celldm(..)` are Bohr.

    Does minimal checking of the cell shape, so it is possible to create
    a nonsense structure if the ibrav is inapproprite for the cell. These
    are derived to be symmetric with the routine for constructing the cell
    from ibrav parameters so directions of some vectors may be unexpected.

    Parameters
    ----------
    cell : np.array
        A 3x3 representation of a unit cell
    ibrav : int
        Bravais-lattice index according to the pw.x designations.

    Returns
    -------
    parameters : dict
        A dictionary with all the necessary `celldm(..)` keys assigned
        necessary values (in units of Bohr). Also includes `ibrav` so it
        can be passed back to `ibrav_to_cell`.

    Raises
    ------
    NotImplementedError
        Only a limited number of ibrav settings can be parsed. An error
        is raised if the ibrav interpretation is not implemented.
    """
    parameters = {'ibrav': ibrav}

    if ibrav == 1:
        parameters['celldm(1)'] = cell[0][0] / units['Bohr']
    elif ibrav in [2, 3, -3]:
        parameters['celldm(1)'] = cell[0][2] * 2 / units['Bohr']
    elif ibrav in [4, 6]:
        parameters['celldm(1)'] = cell[0][0] / units['Bohr']
        parameters['celldm(3)'] = cell[2][2] / cell[0][0]
    elif ibrav in [5, -5]:
        # Manually derive
        a = np.linalg.norm(cell[0])
        cosab = np.dot(cell[0], cell[1]) / (a ** 2)
        parameters['celldm(1)'] = a / units['Bohr']
        parameters['celldm(4)'] = cosab
    elif ibrav == 7:
        parameters['celldm(1)'] = cell[0][0] * 2 / units['Bohr']
        parameters['celldm(3)'] = cell[2][2] / cell[0][0]
    elif ibrav == 8:
        parameters['celldm(1)'] = cell[0][0] / units['Bohr']
        parameters['celldm(2)'] = cell[1][1] / cell[0][0]
        parameters['celldm(3)'] = cell[2][2] / cell[0][0]
    elif ibrav in [9, -9]:
        parameters['celldm(1)'] = cell[0][0] * 2 / units['Bohr']
        parameters['celldm(2)'] = cell[1][1] / cell[0][0]
        parameters['celldm(3)'] = cell[2][2] * 2 / cell[0][0]
    elif ibrav in [10, 11]:
        parameters['celldm(1)'] = cell[0][0] * 2 / units['Bohr']
        parameters['celldm(2)'] = cell[1][1] / cell[0][0]
        parameters['celldm(3)'] = cell[2][2] / cell[0][0]
    elif ibrav == 12:
        # cos^2 + sin^2
        b = (cell[1][0]**2 + cell[1][1]**2)**0.5
        parameters['celldm(1)'] = cell[0][0] / units['Bohr']
        parameters['celldm(2)'] = b / cell[0][0]
        parameters['celldm(3)'] = cell[2][2] / cell[0][0]
        parameters['celldm(4)'] = cell[1][0] / b
    elif ibrav == -12:
        # cos^2 + sin^2
        c = (cell[2][0]**2 + cell[2][2]**2)**0.5
        parameters['celldm(1)'] = cell[0][0] / units['Bohr']
        parameters['celldm(2)'] = cell[1][1] / cell[0][0]
        parameters['celldm(3)'] = c / cell[0][0]
        parameters['celldm(4)'] = cell[2][0] / c
    elif ibrav == 13:
        b = (cell[1][0]**2 + cell[1][1]**2)**0.5
        parameters['celldm(1)'] = cell[0][0] * 2 / units['Bohr']
        parameters['celldm(2)'] = b / (cell[0][0] * 2)
        parameters['celldm(3)'] = cell[2][2] / cell[0][0]
        parameters['celldm(4)'] = cell[1][0] / b
    elif ibrav == 14:
        # Manually derive
        a, b, c = np.linalg.norm(cell, axis=1)
        cosbc = np.dot(cell[1], cell[2]) / (b * c)
        cosac = np.dot(cell[0], cell[2]) / (a * c)
        cosab = np.dot(cell[0], cell[1]) / (a * b)
        parameters['celldm(1)'] = a / units['Bohr']
        parameters['celldm(2)'] = b / a
        parameters['celldm(3)'] = c / a
        parameters['celldm(4)'] = cosbc
        parameters['celldm(5)'] = cosac
        parameters['celldm(6)'] = cosab
    else:
        raise NotImplementedError('ibrav = {0} is not implemented'
                                  ''.format(ibrav))

    return parameters


def kspacing_to_grid(atoms, spacing, calculated_spacing=None):
    """
    Calculate the kpoint mesh that is equivalent to the given spacing
    in reciprocal space (units Angstrom^-1). The number of kpoints is each
    dimension is rounded up (compatible with CASTEP).

    Parameters
    ----------
    atoms: ase.Atoms
        A structure that can have get_reciprocal_cell called on it.
    spacing: float
        Minimum K-Point spacing in $A^{-1}$.
    calculated_spacing : list
        If a three item list (or similar mutable sequence) is given the
        members will be replaced with the actual calculated spacing in
        $A^{-1}$.

    Returns
    -------
    kpoint_grid : [int, int, int]
        MP grid specification to give the required spacing.

    """
    # No factor of 2pi in ase, everything in A^-1
    # reciprocal dimensions
    r_x, r_y, r_z = np.linalg.norm(atoms.get_reciprocal_cell(), axis=1)

    kpoint_grid = [int(r_x / spacing) + 1,
                   int(r_y / spacing) + 1,
                   int(r_z / spacing) + 1]

    if calculated_spacing is not None:
        calculated_spacing[:] = [r_x / kpoint_grid[0],
                                 r_y / kpoint_grid[1],
                                 r_z / kpoint_grid[2]]

    return kpoint_grid


def construct_kpoints_card(atoms, kpts=None, kspacing=None, koffset=(0, 0, 0)):
    out = []
    kpts_via_parameters = atoms.calc.parameters.get('kpts', None)
    if kpts is None and kpts_via_parameters is not None:
        kpts = kpts_via_parameters

    if kspacing is not None:
        kgrid = kspacing_to_grid(atoms, kspacing)
    elif kpts is not None:
        if isinstance(kpts, dict) and 'path' not in kpts:
            kgrid, shift = kpts2sizeandoffsets(atoms=atoms, **kpts)
            koffset = []
            for i, x in enumerate(shift):
                assert x == 0 or abs(x * kgrid[i] - 0.5) < 1e-14
                koffset.append(0 if x == 0 else 1)
        else:
            kgrid = kpts
    else:
        kgrid = "gamma"

    # True and False work here and will get converted by ':d' format
    if isinstance(koffset, int):
        koffset = (koffset, ) * 3

    # BandPath object or bandpath-as-dictionary:
    if isinstance(kgrid, dict) or hasattr(kgrid, 'kpts'):
        out.append('K_POINTS crystal_b\n')
        assert hasattr(kgrid, 'path') or 'path' in kgrid
        kgrid = kpts2ndarray(kgrid, atoms=atoms)
        out.append('%s\n' % len(kgrid))
        for k in kgrid:
            out.append('{k[0]:.14f} {k[1]:.14f} {k[2]:.14f} 0\n'.format(k=k))
        out.append('\n')
    elif isinstance(kgrid, str) and (kgrid == "gamma"):
        out.append('K_POINTS gamma\n')
        out.append('\n')
    else:
        out.append('K_POINTS automatic\n')
        out.append('{0[0]} {0[1]} {0[2]}  {1[0]:d} {1[1]:d} {1[2]:d}\n'
                   ''.format(kgrid, koffset))
        out.append('\n')
    return out


def get_constraint(constraint_idx):
    """
    Map constraints from QE input/output to FixAtoms or FixCartesian constraint
    """
    if not np.any(constraint_idx):
        return None

    a = [a for a, c in enumerate(constraint_idx) if np.all(c is not None)]
    mask = [[(ic + 1) % 2 for ic in c] for c in constraint_idx
            if np.all(c is not None)]

    if np.all(np.array(mask)) == 1:
        constraint = FixAtoms(a)
    else:
        constraint = FixCartesian(a, mask)
    return constraint


def time_to_float(time_str):
    if 'd' in time_str:
        days = time_str.split('d')[0]
        time_str = time_str.split('d')[1].strip()
    else:
        days = 0
    if 'h' in time_str:
        hours, rem = time_str.split('h')
    if 'h' in time_str:
        hours, rem = time_str.split('h')
    else:
        hours, rem = 0, time_str
    if 'm' in rem:
        minutes, rem = rem.split('m')
    else:
        minutes = 0
    if 's' in rem:
        seconds = rem.rstrip('s')
    else:
        seconds = 0
    return (float(days) * 1440 + float(hours) * 60 + float(minutes)) * 60 + float(seconds)


def safe_float(string):
    # Equivalent to float(string), but if string = '*******' it returns np.nan
    if all([c == '*' for c in string]):
        return np.nan
    else:
        return float(string)


def safe_string_to_list_of_floats(string):
    '''
    Converts a string to a list of floats, but if string = 'number*******anothernumber' it returns
    [number, np.nan, anothernumber]. It also keeps its eye out for cases like numberanothernumber
    and number-negativenumber
    '''

    out = []

    string.replace('-', ' -')
    for word in string.split():
        if '*' in word:
            # First, reduce each sequence of '*'s to a single '*'
            while '**' in word:
                word = word.replace('**', '*')

            # Then pad each '*' with spaces, making sure we don't add spaces at the start/end
            word = word.replace('*', ' * ').strip()

            # Finally, convert to float
            out += [safe_float(x) for x in word.split()]
        else:
            if word.count('.') > 1:
                # For cases where a number almost overflows we might get a single "word" with multiple
                # decimal points...
                out += [np.nan for _ in range(word.count('.'))]
            else:
                out.append(safe_float(word))
    return out


def write_espresso_in(fd, atoms, input_data=None, pseudopotentials=None,
                      kspacing=None, kpts=None, koffset=(0, 0, 0),
                      crystal_coordinates=None, local_construct_namelist=None,
                      include_kpoints=True, gamma_only=False, **kwargs):
    """
    Create an input file for a generic Quantum ESPRESSO calculator

    Use set_initial_magnetic_moments to turn on spin, if ispin is set to 2
    with no magnetic moments, they will all be set to 0.0. Magnetic moments
    will be converted to the QE units (fraction of valence electrons) using
    any pseudopotential files found, or a best guess for the number of
    valence electrons.

    Units are not converted for any other input data, so use Quantum ESPRESSO
    units (Usually Ry or atomic units).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    Implemented features:

    - Conversion of :class:`ase.constraints.FixAtoms` and
                    :class:`ase.constraints.FixCartesian`.
    - `starting_magnetization` derived from the `mgmoms` and pseudopotentials
      (searches default paths for pseudo files.)
    - Automatic assignment of options to their correct sections.
    - Interpretation of ibrav (cell must exactly match the vectors defined
      in the QE docs).

    Not implemented:

    - Lists of k-points
    - Other constraints
    - Hubbard parameters
    - Validation of the argument types for input
    - Validation of required options
    - Reorientation for ibrav settings
    - Noncollinear magnetism

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.
    input_data: dict
        A flat or nested dictionary with input parameters for pw.x
    pseudopotentials: dict
        A filename for each atomic species, e.g.
        {'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}.
        A dummy name will be used if none are given.
    kspacing: float
        Generate a grid of k-points with this as the minimum distance,
        in A^-1 between them in reciprocal space. If set to None, kpts
        will be used instead.
    kpts: (int, int, int) or dict
        If kpts is a tuple (or list) of 3 integers, it is interpreted
        as the dimensions of a Monkhorst-Pack grid.
        If ``kpts`` is set to ``None``, only the Γ-point will be included
        and QE will use routines optimized for Γ-point-only calculations.
        Compared to Γ-point-only calculations without this optimization
        (i.e. with ``kpts=(1, 1, 1)``), the memory and CPU requirements
        are typically reduced by half.
        If kpts is a dict, it will either be interpreted as a path
        in the Brillouin zone (*) if it contains the 'path' keyword,
        otherwise it is converted to a Monkhorst-Pack grid (**).
        (*) see ase.dft.kpoints.bandpath
        (**) see ase.calculators.calculator.kpts2sizeandoffsets
    koffset: (int, int, int)
        Offset of kpoints in each direction. Must be 0 (no offset) or
        1 (half grid offset). Setting to True is equivalent to (1, 1, 1).
    crystal_coordinates: bool
        Whether the atomic positions should be written to the QE input file in
        absolute (False) or relative (crystal) coordinates (True). If no value
        is provided, this keyword will take the value of atoms.pbc

    """

    # Default for crystal_coordinates
    if crystal_coordinates is None:
        crystal_coordinates = all(atoms.pbc)

    # Convert to a namelist to make working with parameters much easier
    # Note that the name ``input_data`` is chosen to prevent clash with
    # ``parameters`` in Calculator objects
    if 'input_data' in atoms.calc.parameters:
        if input_data is None:
            input_data = atoms.calc.parameters['input_data']
        elif input_data != atoms.calc.parameters['input_data']:
            warnings.warn("write_espresso_in(...) is ignoring "
                          "atom.calc.parameters['input_data'] in favor of "
                          "the input_data provided to it as an argument")

    input_parameters = local_construct_namelist(input_data, **kwargs)

    # Convert ase constraints to QE constraints
    # Nx3 array of force multipliers matches what QE uses
    # Do this early so it is available when constructing the atoms card
    constraint_mask = np.ones((len(atoms), 3), dtype='int')
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            constraint_mask[constraint.index] = 0
        elif isinstance(constraint, FixCartesian):
            constraint_mask[constraint.a] = constraint.mask
        else:
            warnings.warn('Ignored unknown constraint {}'.format(constraint))

    # Deal with pseudopotentials
    # Look in all possible locations for the pseudos and try to figure
    # out the number of valence electrons
    pseudo_dirs = []
    if 'pseudo_dir' in input_parameters['control']:
        pseudo_dirs.append(input_parameters['control']['pseudo_dir'])
    if 'ESPRESSO_PSEUDO' in os.environ:
        pseudo_dirs.append(os.environ['ESPRESSO_PSEUDO'])
    pseudo_dirs.append(os.path.expanduser('~/espresso/pseudo/'))

    # Species info holds the information on the pseudopotential and
    # associated for each element
    if pseudopotentials is None:
        pseudopotentials = atoms.calc.parameters.get('pseudopotentials', {})
    species_info = {}

    species = atoms.get_chemical_symbols()
    if len(set(atoms.get_tags())) > 1:
        labels = [s + str(t) if t > 0 else s for s, t in zip(species, atoms.get_tags())]
    else:
        labels = species

    for label, specie in set(zip(labels, species)):
        pseudo = pseudopotentials.get(label, '{}_dummy.UPF'.format(specie))
        for pseudo_dir in pseudo_dirs:
            if os.path.exists(os.path.join(pseudo_dir, pseudo)):
                valence = grep_valence(os.path.join(pseudo_dir, pseudo))
                break
        else:  # not found in a file
            valence = SSSP_VALENCE[atomic_numbers[specie]]

        species_info[label] = {'pseudo': pseudo,
                               'valence': valence}

    # Convert atoms into species.
    # Each different magnetic moment needs to be a separate type even with
    # the same pseudopotential (e.g. an up and a down for AFM).
    # if any magmom are > 0 or nspin == 2 then use species labels.
    # Rememeber: magnetisation uses 1 based indexes
    atomic_species = OrderedDict()
    atomic_species_str = []
    atomic_positions_str = []

    nspin = input_parameters['system'].get('nspin', 1)  # 1 is the default
    if any(atoms.get_initial_magnetic_moments()):
        if nspin == 1:
            # Force spin on
            input_parameters['system']['nspin'] = 2
            nspin = 2

    if nspin == 2:
        # Spin on
        for atom, label, magmom in zip(atoms, labels, atoms.get_initial_magnetic_moments()):
            if (label, magmom) not in atomic_species:
                # spin as fraction of valence
                fspin = float(magmom) / species_info[label]['valence']
                # Index in the atomic species list
                sidx = len(atomic_species) + 1
                # Index for that atom type; no index for first one
                tidx = sum(atom.symbol == x[0] for x in atomic_species) or ' '
                atomic_species[(label, magmom)] = (sidx, tidx)
                # Add magnetization to the input file
                mag_str = 'starting_magnetization({0})'.format(sidx)
                input_parameters['system'][mag_str] = fspin
                atomic_species_str.append(
                    '{species}{tidx} {mass} {pseudo}\n'.format(
                        species=label, tidx=tidx, mass=atom.mass,
                        pseudo=species_info[label]['pseudo']))
            # lookup tidx to append to name
            sidx, tidx = atomic_species[(label, magmom)]

            # only inclued mask if something is fixed
            if not all(constraint_mask[atom.index]):
                mask = ' {mask[0]} {mask[1]} {mask[2]}'.format(
                    mask=constraint_mask[atom.index])
            else:
                mask = ''

            # construct line for atomic positions
            if crystal_coordinates:
                coords = [atom.a, atom.b, atom.c]
            else:
                coords = atom.position
            atomic_positions_str.append(
                f'{label}{tidx} '
                f'{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}'
                f'{mask}\n')

    else:
        # Wipe out tot_magnetisation
        if input_parameters['system'].pop('tot_magnetization', 0) != 0:
            raise ValueError('tot_magnetization != 0 and nspin == 1 are incompatible')

        # Do nothing about magnetisation
        for atom, label in zip(atoms, labels):
            if label not in atomic_species:
                atomic_species[label] = True  # just a placeholder
                atomic_species_str.append(
                    '{species} {mass} {pseudo}\n'.format(
                        species=label, mass=atom.mass,
                        pseudo=species_info[label]['pseudo']))

            # only inclued mask if something is fixed
            if not all(constraint_mask[atom.index]):
                mask = ' {mask[0]} {mask[1]} {mask[2]}'.format(
                    mask=constraint_mask[atom.index])
            else:
                mask = ''

            if crystal_coordinates:
                coords = [atom.a, atom.b, atom.c]
            else:
                coords = atom.position
            atomic_positions_str.append(f'{label} {coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f} {mask}\n')

        if input_parameters['system'].pop('tot_magnetization', 0) != 0:
            raise ValueError('tot_magnetization cannot be non-zero when nspin = 1')

    # Add computed parameters
    # different magnetisms means different types
    input_parameters['system']['ntyp'] = len(atomic_species)
    input_parameters['system']['nat'] = len(atoms)

    # Use cell as given or fit to a specific ibrav
    if 'ibrav' in input_parameters['system']:
        ibrav = input_parameters['system']['ibrav']
        if ibrav != 0:
            celldm = cell_to_ibrav(atoms.cell, ibrav)
            regen_cell = ibrav_to_cell(celldm)[1]
            if not np.allclose(atoms.cell, regen_cell):
                warnings.warn('Input cell does not match requested ibrav'
                              '{} != {}'.format(regen_cell, atoms.cell))
            input_parameters['system'].update(celldm)
    else:
        # Just use standard cell block
        input_parameters['system']['ibrav'] = 0

    # Construct input file into this
    flines = []

    # Assume sections are ordered (taken care of in namelist construction)
    # and that repr converts to a QE readable representation (except bools)
    for section in input_parameters:
        flines.append('&{0}\n'.format(section.upper()))
        for key, value in input_parameters[section].items():
            if value is True:
                flines.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                flines.append('   {0:16} = .false.\n'.format(key))
            elif value is not None:
                # repr format to get quotes around strings
                flines.append('   {0:16} = {1!r:}\n'.format(key, value))
        flines.append('/\n')  # terminate section
    flines.append('\n')

    flines.append('ATOMIC_SPECIES\n')
    flines.extend(atomic_species_str)
    flines.append('\n')

    # KPOINTS - add a MP grid as required
    if include_kpoints:
        if gamma_only:
            flines += construct_kpoints_card(atoms, None, None)
        else:
            flines += construct_kpoints_card(atoms, kpts, kspacing, koffset)

    # CELL block, if required
    if input_parameters['SYSTEM']['ibrav'] == 0:
        flines.append('CELL_PARAMETERS angstrom\n')
        flines.append('{cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
                      '{cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
                      '{cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}\n'
                      ''.format(cell=atoms.cell))
        flines.append('\n')

    # Positions - already constructed, but must appear after namelist
    if crystal_coordinates:
        flines.append('ATOMIC_POSITIONS crystal\n')
    else:
        flines.append('ATOMIC_POSITIONS angstrom\n')
    flines.extend(atomic_positions_str)
    flines.append('\n')

    # DONE!
    fd.write(''.join(flines))
