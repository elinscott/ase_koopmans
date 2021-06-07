""" Utilities for Quantum ESPRESSO file I/O

Written by Edward Linscott 2020-21

"""

import warnings
import numpy as np
import operator as op
from collections import OrderedDict
from ase.units import create_units
from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
from ase.constraints import FixAtoms, FixCartesian
from ase.data import chemical_symbols

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


def get_kpoints(card_lines):
    # Find kpts, koffset from card_lines
    for i, line in enumerate(card_lines):
        if line.startswith('K_POINTS'):
            mode = line.split()[-1]
            if mode.lower() == 'gamma':
                return [1, 1, 1], [0, 0, 0]
            elif mode.lower() == 'automatic':
                splitline = card_lines[i + 1].strip().split()
                kpts = [int(x) for x in splitline[:3]]
                koffset = [int(x) for x in splitline[3:]]
                return kpts, koffset
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
                if alat is not None:
                    raise ValueError('Lattice parameters given in '
                                     '&SYSTEM celldm/A and CELL_PARAMETERS '
                                     'angstrom')
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


def label_to_symbol(label):
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

    Raises
    ------
    KeyError
        Couldn't find an appropriate species.

    Notes
    -----
        It's impossible to tell whether e.g. He is helium
        or hydrogen labelled 'e'.
    """

    # possibly a two character species
    # ase Atoms need proper case of chemical symbols.
    if len(label) >= 2:
        test_symbol = label[0].upper() + label[1].lower()
        if test_symbol in chemical_symbols:
            return test_symbol
    # finally try with one character
    test_symbol = label[0].upper()
    if test_symbol in chemical_symbols:
        return test_symbol
    else:
        raise KeyError('Could not parse species from label {0}.'
                       ''.format(label))


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
            if key in parameters.get(section, {}):
                sec_list[key] = parameters[section].pop(key)
            if key in parameters:
                sec_list[key] = parameters.pop(key)
            if key in kwargs:
                sec_list[key] = kwargs.pop(key)

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
    # Converts a string to a list of floats, but if string = 'number*******anothernumber' it returns [number, np.nan, anothernumber]
    out = []
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
            out.append(safe_float(word))
    return out
