import os
import re
from ase.utils import reader


"""ASE-interface to Octopus.

Ask Hjorth Larsen <asklarsen@gmail.com>
Carlos de Armas

http://tddft.org/programs/octopus/
"""

import numpy as np

from ase import Atoms
from ase.data import atomic_numbers
from ase.io import read
from ase.units import Bohr, Angstrom, Hartree, eV, Debye
from ase.calculators.calculator import kpts2ndarray


# Representation of parameters from highest to lowest level of abstraction:
#
#  * Atoms object plus reduced kwargs that specify info not stored in the Atoms
#  * full dictionary of kwargs incorporating info contained in the Atoms
#  * series of assignments (names, values).  May contain duplicates.
#  * text in Octopus input file format


# Octopus variable types and specification in Python:
#
#  Type     Examples                    Equivalent in Python:
# -----------------------------------------------------------------------------
#  flag     wfs + density               'wfs + density'
#  float    2.7, 2.7 + pi^2             2.7, 2.7 + np.pi**2
#  integer  42, rmmdiis                 42, 'rmmdiis'
#  logical  true, false, yes, no, ...   True, False, 1, 0, 'yes', 'no', ...
#  string   "stdout.txt"                '"stdout.txt"' (apologies for ugliness)
#
#  block    %Coordinates                List of lists:
#            'H' | 0 | 0 | 0              coordinates=[["'H'", 0, 0, 0],
#            'O' | 0 | 0 | 1                           ["'O'", 0, 0, 1]]
#           %                             (elemements are sent through repr())

# Rules for input parameters
# --------------------------
#
# We make the following conversions:
#     dict of keyword arguments + Atoms object -> Octopus input file
# and:
#     Octopus input file -> Atoms object + dict of keyword arguments
# Below we detail some conventions and compatibility issues.
#
# 1) ASE always passes some parameters by default (always write
#    forces, etc.).  They can be overridden by the user, but the
#    resulting behaviour is undefined.
#
# 2) Atoms object is used to establish some parameters: Coordinates,
#    Lsize, etc.  All those parameters can be overridden by passing
#    them directly as keyword arguments.  Parameters that were taken
#    from the Atoms object are always marked with the comment "# ASE
#    auto" in the input file.  This is used to distinguish variables
#    that are overridden from variables that simply came from the
#    atoms object when restarting.
#
# 3) Some variables do not interact nicely between ASE and Octopus,
#    such as SubSystemCoordinates which may involve rotations.  There
#    may be many such variables that we have not identified, but at
#    least the known ones will cause a suppressable
#    OctopusKeywordError.  (This third rule has not been implemented
#    as of this moment.)
#
# 4) OctopusKeywordError is raised from Python for keywords that are
#    not valid according to oct-help.


def read_eigenvalues_file(fd):
    unit = None

    for line in fd:
        m = re.match(r'Eigenvalues\s*\[(.+?)\]', line)
        if m is not None:
            unit = m.group(1)
            break
    line = next(fd)
    assert line.strip().startswith('#st'), line

    #fermilevel = None
    kpts = []
    eigs = []
    occs = []

    for line in fd:
        m = re.match(r'#k.*?\(\s*(.+?),\s*(.+?),\s*(.+?)\)', line)
        if m:
            k = m.group(1, 2, 3)
            kpts.append(np.array(k, float))
            eigs.append({})
            occs.append({})
        else:
            m = re.match(r'\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)', line)
            if m is None:
                m = re.match(r'Fermi energy\s*=\s*(\S+)\s*', line)
                assert m is not None
                # We can also return the fermilevel but so far we just read
                # it from the static/info instead.
                #fermilevel = float(m.group(1))
            else:
                spin, eig, occ = m.group(1, 2, 3)

                if not eigs:
                    # Only initialized if kpoint header was written
                    eigs.append({})
                    occs.append({})

                eigs[-1].setdefault(spin, []).append(float(eig))
                occs[-1].setdefault(spin, []).append(float(occ))

    nkpts = len(kpts)
    nspins = len(eigs[0])
    nbands = len(eigs[0][spin])

    kptsarr = np.array(kpts, float)
    eigsarr = np.empty((nkpts, nspins, nbands))
    occsarr = np.empty((nkpts, nspins, nbands))

    arrs = [eigsarr, occsarr]

    for arr in arrs:
        arr.fill(np.nan)

    for k in range(nkpts):
        for arr, lst in [(eigsarr, eigs), (occsarr, occs)]:
            arr[k, :, :] = [lst[k][sp] for sp
                            in (['--'] if nspins == 1 else ['up', 'dn'])]

    for arr in arrs:
        assert not np.isnan(arr).any()

    eigsarr *= {'H': Hartree, 'eV': eV}[unit]
    return kptsarr, eigsarr, occsarr


def process_special_kwargs(atoms, kwargs):
    kwargs = kwargs.copy()
    kpts = kwargs.pop('kpts', None)
    if kpts is not None:
        for kw in ['kpoints', 'reducedkpoints', 'kpointsgrid']:
            if kw in kwargs:
                raise ValueError('k-points specified multiple times')

        kptsarray = kpts2ndarray(kpts, atoms)
        nkpts = len(kptsarray)
        fullarray = np.empty((nkpts, 4))
        fullarray[:, 0] = 1.0 / nkpts  # weights
        fullarray[:, 1:4] = kptsarray
        kwargs['kpointsreduced'] = fullarray.tolist()

    # TODO xc=LDA/PBE etc.

    # The idea is to get rid of the special keywords, since the rest
    # will be passed to Octopus
    # XXX do a better check of this
    from ase.calculators.octopus import Octopus
    for kw in Octopus.special_ase_keywords:
        assert kw not in kwargs, kw
    return kwargs


class OctopusKeywordError(ValueError):
    pass  # Unhandled keywords


class OctopusParseError(Exception):
    pass  # Cannot parse input file


class OctopusIOError(IOError):
    pass  # Cannot find output files


def unpad(pbc, arr):
    # Return non-padded array from padded array.
    # This means removing the last element along all periodic directions.
    if pbc[0]:
        assert np.all(arr[0, :, :] == arr[-1, :, :])
        arr = arr[0:-1, :, :]
    if pbc[1]:
        assert np.all(arr[:, 0, :] == arr[:, -1, :])
        arr = arr[:, 0:-1, :]
    if pbc[2]:
        assert np.all(arr[:, :, 0] == arr[:, :, -1])
        arr = arr[:, :, 0:-1]
    return np.ascontiguousarray(arr)


def unpad_smarter(pbc, arr):
    # 'Smarter' but less easy to understand version of the above.
    # (untested I think)
    slices = []
    for c, is_periodic in enumerate(pbc):
        if is_periodic:
            left = np.take(arr, [0], axis=c)
            right = np.take(arr, [-1], axis=c)
            assert np.all(left == right)
            slices.append(slice(0, -1))
        else:
            slices.append(slice(None))
    return np.ascontiguousarray(arr[slices])


# Parse value as written in input file *or* something that one would be
# passing to the ASE interface, i.e., this might already be a boolean
def octbool2bool(value):
    value = value.lower()
    if isinstance(value, int):
        return bool(value)
    if value in ['true', 't', 'yes', '1']:
        return True
    elif value in ['no', 'f', 'false', '0']:
        return False
    else:
        raise ValueError('Failed to interpret "%s" as a boolean.' % value)


def list2block(name, rows):
    """Construct 'block' of Octopus input.

    convert a list of rows to a string with the format x | x | ....
    for the octopus input file"""
    lines = []
    lines.append('%' + name)
    for row in rows:
        lines.append(' ' + ' | '.join(str(obj) for obj in row))
    lines.append('%')
    return lines


def normalize_keywords(kwargs):
    """Reduce keywords to unambiguous form (lowercase)."""
    newkwargs = {}
    for arg, value in kwargs.items():
        lkey = arg.lower()
        newkwargs[lkey] = value
    return newkwargs


def input_line_iter(lines):
    """Convenient iterator for parsing input files 'cleanly'.

    Discards comments etc."""
    for line in lines:
        line = line.split('#')[0].strip()
        if not line or line.isspace():
            continue
        line = line.strip()
        yield line


def block2list(namespace, lines, header=None):
    """Parse lines of block and return list of lists of strings."""
    lines = iter(lines)
    block = []
    if header is None:
        header = next(lines)
    assert header.startswith('%'), header
    name = header[1:].strip().lower()
    for line in lines:
        if line.startswith('%'):  # Could also say line == '%' most likely.
            break
        tokens = [namespace.evaluate(token)
                  for token in line.strip().split('|')]
        # XXX will fail for string literals containing '|'
        block.append(tokens)
    return name, block


class OctNamespace:
    def __init__(self):
        self.names = {}
        self.consts = {'pi': np.pi,
                       'angstrom': 1. / Bohr,
                       'ev': 1. / Hartree,
                       'yes': True,
                       'no': False,
                       't': True,
                       'f': False,
                       'i': 1j,  # This will probably cause trouble
                       'true': True,
                       'false': False}

    def evaluate(self, value):
        value = value.strip()

        for char in '"', "'":  # String literal
            if value.startswith(char):
                assert value.endswith(char)
                return value

        value = value.lower()

        if value in self.consts:  # boolean or other constant
            return self.consts[value]

        if value in self.names:  # existing variable
            return self.names[value]

        try:  # literal integer
            v = int(value)
        except ValueError:
            pass
        else:
            if v == float(v):
                return v

        try:  # literal float
            return float(value)
        except ValueError:
            pass

        if ('*' in value or '/' in value
            and not any(char in value for char in '()+')):
            floatvalue = 1.0
            op = '*'
            for token in re.split(r'([\*/])', value):
                if token in '*/':
                    op = token
                    continue

                v = self.evaluate(token)

                try:
                    v = float(v)
                except TypeError:
                    try:
                        v = complex(v)
                    except ValueError:
                        break
                except ValueError:
                    break  # Cannot evaluate expression
                else:
                    if op == '*':
                        floatvalue *= v
                    else:
                        assert op == '/', op
                        floatvalue /= v
            else:  # Loop completed successfully
                return floatvalue
        return value  # unknown name, or complex arithmetic expression

    def add(self, name, value):
        value = self.evaluate(value)
        self.names[name.lower().strip()] = value


def parse_input_file(fd):
    namespace = OctNamespace()
    lines = input_line_iter(fd)
    blocks = {}
    while True:
        try:
            line = next(lines)
        except StopIteration:
            break
        else:
            if line.startswith('%'):
                name, value = block2list(namespace, lines, header=line)
                blocks[name] = value
            else:
                tokens = line.split('=', 1)
                assert len(tokens) == 2, tokens
                name, value = tokens
                namespace.add(name, value)

    namespace.names.update(blocks)
    return namespace.names


def kwargs2cell(kwargs):
    # kwargs -> cell + remaining kwargs
    # cell will be None if not ASE-compatible.
    #
    # Returns numbers verbatim; caller must convert units.
    kwargs = normalize_keywords(kwargs)

    if boxshape_is_ase_compatible(kwargs):
        kwargs.pop('boxshape', None)
        if 'lsize' in kwargs:
            Lsize = kwargs.pop('lsize')
            if not isinstance(Lsize, list):
                Lsize = [[Lsize] * 3]
            assert len(Lsize) == 1
            cell = np.array([2 * float(l) for l in Lsize[0]])
        elif 'latticeparameters' in kwargs:
            # Eval latparam and latvec
            latparam = np.array(kwargs.pop('latticeparameters'), float).T
            cell = np.array(kwargs.pop('latticevectors', np.eye(3)), float)
            for a, vec in zip(latparam, cell):
                vec *= a
            assert cell.shape == (3, 3)
    else:
        cell = None
    return cell, kwargs


def boxshape_is_ase_compatible(kwargs):
    pdims = int(kwargs.get('periodicdimensions', 0))
    default_boxshape = 'parallelepiped' if pdims > 0 else 'minimum'
    boxshape = kwargs.get('boxshape', default_boxshape).lower()
    # XXX add support for experimental keyword 'latticevectors'
    return boxshape == 'parallelepiped'


def kwargs2atoms(kwargs, directory=None):
    """Extract atoms object from keywords and return remaining keywords.

    Some keyword arguments may refer to files.  The directory keyword
    may be necessary to resolve the paths correctly, and is used for
    example when running 'ase gui somedir/inp'."""
    kwargs = normalize_keywords(kwargs)

    # Only input units accepted nowadays are 'atomic'.
    # But if we are loading an old file, and it specifies something else,
    # we can be sure that the user wanted that back then.

    coord_keywords = ['coordinates',
                      'xyzcoordinates',
                      'pdbcoordinates',
                      'reducedcoordinates',
                      'xsfcoordinates',
                      'xsfcoordinatesanimstep']

    nkeywords = 0
    for keyword in coord_keywords:
        if keyword in kwargs:
            nkeywords += 1
    if nkeywords == 0:
        raise OctopusParseError('No coordinates')
    elif nkeywords > 1:
        raise OctopusParseError('Multiple coordinate specifications present.  '
                                'This may be okay in Octopus, but we do not '
                                'implement it.')

    def get_positions_from_block(keyword):
        # %Coordinates or %ReducedCoordinates -> atomic numbers, positions.
        block = kwargs.pop(keyword)
        positions = []
        numbers = []
        tags = []
        types = {}
        for row in block:
            assert len(row) in [ndims + 1, ndims + 2]
            row = row[:ndims + 1]
            sym = row[0]
            assert sym.startswith('"') or sym.startswith("'")
            assert sym[0] == sym[-1]
            sym = sym[1:-1]
            pos0 = np.zeros(3)
            ndim = int(kwargs.get('dimensions', 3))
            pos0[:ndim] = [float(element) for element in row[1:]]
            number = atomic_numbers.get(sym)  # Use 0 ~ 'X' for unknown?
            tag = 0
            if number is None:
                if sym not in types:
                    tag = len(types) + 1
                    types[sym] = tag
                number = 0
                tag = types[sym]
            tags.append(tag)
            numbers.append(number)
            positions.append(pos0)
        positions = np.array(positions)
        tags = np.array(tags, int)
        if types:
            ase_types = {}
            for sym, tag in types.items():
                ase_types[('X', tag)] = sym
            info = {'types': ase_types}  # 'info' dict for Atoms object
        else:
            tags = None
            info = None
        return numbers, positions, tags, info

    def read_atoms_from_file(fname, fmt):
        assert fname.startswith('"') or fname.startswith("'")
        assert fname[0] == fname[-1]
        fname = fname[1:-1]
        if directory is not None:
            fname = os.path.join(directory, fname)
        # XXX test xyz, pbd and xsf
        if fmt == 'xsf' and 'xsfcoordinatesanimstep' in kwargs:
            anim_step = kwargs.pop('xsfcoordinatesanimstep')
            theslice = slice(anim_step, anim_step + 1, 1)
            # XXX test animstep
        else:
            theslice = slice(None, None, 1)
        images = read(fname, theslice, fmt)
        if len(images) != 1:
            raise OctopusParseError('Expected only one image.  Don\'t know '
                                    'what to do with %d images.' % len(images))
        return images[0]

    # We will attempt to extract cell and pbc from kwargs if 'lacking'.
    # But they might have been left unspecified on purpose.
    #
    # We need to keep track of these two variables "externally"
    # because the Atoms object assigns values when they are not given.
    cell = None
    pbc = None
    adjust_positions_by_half_cell = False

    atoms = None
    xsfcoords = kwargs.pop('xsfcoordinates', None)
    if xsfcoords is not None:
        atoms = read_atoms_from_file(xsfcoords, 'xsf')
        atoms.positions *= Bohr
        atoms.cell *= Bohr
        # As it turns out, non-periodic xsf is not supported by octopus.
        # Also, it only supports fully periodic or fully non-periodic....
        # So the only thing that we can test here is 3D fully periodic.
        if sum(atoms.pbc) != 3:
            raise NotImplementedError('XSF not fully periodic with Octopus')
        cell = atoms.cell
        pbc = atoms.pbc
        # Position adjustment doesn't actually matter but this should work
        # most 'nicely':
        adjust_positions_by_half_cell = False
    xyzcoords = kwargs.pop('xyzcoordinates', None)
    if xyzcoords is not None:
        atoms = read_atoms_from_file(xyzcoords, 'xyz')
        atoms.positions *= Bohr
        adjust_positions_by_half_cell = True
    pdbcoords = kwargs.pop('pdbcoordinates', None)
    if pdbcoords is not None:
        atoms = read_atoms_from_file(pdbcoords, 'pdb')
        pbc = atoms.pbc
        adjust_positions_by_half_cell = True
        # Due to an error in ASE pdb, we can only test the nonperiodic case.
        # atoms.cell *= Bohr # XXX cell?  Not in nonperiodic case...
        atoms.positions *= Bohr
        if sum(atoms.pbc) != 0:
            raise NotImplementedError('Periodic pdb not supported by ASE.')

    if cell is None:
        # cell could not be established from the file, so we set it on the
        # Atoms now if possible:
        cell, kwargs = kwargs2cell(kwargs)
        if cell is not None:
            cell *= Bohr
        if cell is not None and atoms is not None:
            atoms.cell = cell
        # In case of boxshape = sphere and similar, we still do not have
        # a cell.

    ndims = int(kwargs.get('dimensions', 3))
    if ndims != 3:
        raise NotImplementedError('Only 3D calculations supported.')

    coords = kwargs.get('coordinates')
    if coords is not None:
        numbers, pos, tags, info = get_positions_from_block('coordinates')
        pos *= Bohr
        adjust_positions_by_half_cell = True
        atoms = Atoms(cell=cell, numbers=numbers, positions=pos,
                      tags=tags, info=info)
    rcoords = kwargs.get('reducedcoordinates')
    if rcoords is not None:
        numbers, spos, tags, info = get_positions_from_block(
            'reducedcoordinates')
        if cell is None:
            raise ValueError('Cannot figure out what the cell is, '
                             'and thus cannot interpret reduced coordinates.')
        atoms = Atoms(cell=cell, numbers=numbers, scaled_positions=spos,
                      tags=tags, info=info)
    if atoms is None:
        raise OctopusParseError('Apparently there are no atoms.')

    # Either we have non-periodic BCs or the atoms object already
    # got its BCs from reading the file.  In the latter case
    # we shall override only if PeriodicDimensions was given specifically:

    if pbc is None:
        pdims = int(kwargs.pop('periodicdimensions', 0))
        pbc = np.zeros(3, dtype=bool)
        pbc[:pdims] = True
        atoms.pbc = pbc

    if (cell is not None and cell.shape == (3,)
        and adjust_positions_by_half_cell):
        nonpbc = (atoms.pbc == 0)
        atoms.positions[:, nonpbc] += np.array(cell)[None, nonpbc] / 2.0

    return atoms, kwargs


def atoms2kwargs(atoms, use_ase_cell):
    kwargs = {}

    positions = atoms.positions / Bohr

    if use_ase_cell:
        cell = atoms.cell / Bohr
        cell_offset = 0.5 * cell.sum(axis=0)
        positions -= cell_offset
        if atoms.cell.orthorhombic:
            Lsize = 0.5 * np.diag(cell)
            kwargs['lsize'] = [[repr(size) for size in Lsize]]
            # ASE uses (0...cell) while Octopus uses -L/2...L/2.
            # Lsize is really cell / 2, and we have to adjust our
            # positions by subtracting Lsize (see construction of the coords
            # block) in non-periodic directions.
        else:
            kwargs['latticevectors'] = cell.tolist()

    types = atoms.info.get('types', {})

    coord_block = []
    for sym, pos, tag in zip(atoms.get_chemical_symbols(),
                             positions, atoms.get_tags()):
        if sym == 'X':
            sym = types.get((sym, tag))
            if sym is None:
                raise ValueError('Cannot represent atom X without tags and '
                                 'species info in atoms.info')
        coord_block.append([repr(sym)] + [repr(x) for x in pos])

    kwargs['coordinates'] = coord_block
    npbc = sum(atoms.pbc)
    for c in range(npbc):
        if not atoms.pbc[c]:
            msg = ('Boundary conditions of Atoms object inconsistent '
                   'with requirements of Octopus.  pbc must be either '
                   '000, 100, 110, or 111.')
            raise ValueError(msg)
    kwargs['periodicdimensions'] = npbc

    # TODO InitialSpins
    #
    # TODO can use maximumiterations + output/outputformat to extract
    # things from restart file into output files without trouble.
    #
    # Velocities etc.?
    return kwargs


def generate_input(atoms, kwargs):
    """Convert atoms and keyword arguments to Octopus input file."""
    _lines = []

    def append(line):
        _lines.append(line)

    def extend(lines):
        _lines.extend(lines)
        append('')

    def setvar(key, var):
        append('%s = %s' % (key, var))

    for kw in ['lsize', 'latticevectors', 'latticeparameters']:
        assert kw not in kwargs

    defaultboxshape = 'parallelepiped' if atoms.pbc.any() else 'minimum'
    boxshape = kwargs.get('boxshape', defaultboxshape).lower()
    use_ase_cell = (boxshape == 'parallelepiped')
    atomskwargs = atoms2kwargs(atoms, use_ase_cell)

    if use_ase_cell:
        if 'lsize' in atomskwargs:
            block = list2block('LSize', atomskwargs['lsize'])
        elif 'latticevectors' in atomskwargs:
            extend(list2block('LatticeParameters', [[1., 1., 1.]]))
            block = list2block('LatticeVectors', atomskwargs['latticevectors'])
        extend(block)

    # Allow override or issue errors?
    pdim = 'periodicdimensions'
    if pdim in kwargs:
        if int(kwargs[pdim]) != int(atomskwargs[pdim]):
            raise ValueError('Cannot reconcile periodicity in input '
                             'with that of Atoms object')
    setvar('periodicdimensions', atomskwargs[pdim])

    # We like to output forces
    if 'output' in kwargs:
        output_string = kwargs.pop('output')
        output_tokens = [token.strip()
                         for token in output_string.lower().split('+')]
    else:
        output_tokens = []

    if 'forces' not in output_tokens:
        output_tokens.append('forces')
    setvar('output', ' + '.join(output_tokens))
    # It is illegal to have output forces without any OutputFormat.
    # Even though the forces are written in the same format no matter
    # OutputFormat.  Thus we have to make one up:

    if 'outputformat' not in kwargs:
        setvar('outputformat', 'xcrysden')

    for key, val in kwargs.items():
        # Most datatypes are straightforward but blocks require some attention.
        if isinstance(val, list):
            append('')
            dict_data = list2block(key, val)
            extend(dict_data)
        else:
            setvar(key, str(val))
    append('')

    coord_block = list2block('Coordinates', atomskwargs['coordinates'])
    extend(coord_block)
    return '\n'.join(_lines)


def read_static_info_stress(fd):
    stress_cv = np.empty((3, 3))

    headers = next(fd)
    assert headers.strip().startswith('T_{ij}')
    for i in range(3):
        line = next(fd)
        tokens = line.split()
        vec = np.array(tokens[1:4]).astype(float)
        stress_cv[i] = vec
    return stress_cv

def read_static_info_kpoints(fd):
    for line in fd:
        if line.startswith('List of k-points'):
            break

    tokens = next(fd).split()
    assert tokens == ['ik', 'k_x', 'k_y', 'k_z', 'Weight']
    bar = next(fd)
    assert bar.startswith('---')

    kpts = []
    weights = []

    for line in fd:
        # Format:        index   kx      ky      kz     weight
        m = re.match(r'\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)', line)
        if m is None:
            break
        kxyz = m.group(1, 2, 3)
        weight = m.group(4)
        kpts.append(kxyz)
        weights.append(weight)

    ibz_k_points = np.array(kpts, float)
    k_point_weights = np.array(weights, float)
    return dict(ibz_k_points=ibz_k_points, k_point_weights=k_point_weights)


def read_static_info_eigenvalues(fd, energy_unit):

    values_sknx = {}

    nbands = 0
    fermilevel = None
    for line in fd:
        line = line.strip()
        if line.startswith('#'):
            continue
        if not line[:1].isdigit():
            m = re.match(r'Fermi energy\s*=\s*(\S+)', line)
            if m is not None:
                fermilevel = float(m.group(1)) * energy_unit
            break

        tokens = line.split()
        nbands = max(nbands, int(tokens[0]))
        energy = float(tokens[2]) * energy_unit
        occupation = float(tokens[3])
        values_sknx.setdefault(tokens[1], []).append((energy, occupation))

    nspins = len(values_sknx)
    if nspins == 1:
        val = [values_sknx['--']]
    else:
        val = [values_sknx['up'], values_sknx['dn']]
    val = np.array(val, float)
    nkpts, remainder = divmod(len(val[0]), nbands)
    assert remainder == 0

    eps_skn = val[:, :, 0].reshape(nspins, nkpts, nbands)
    occ_skn = val[:, :, 1].reshape(nspins, nkpts, nbands)
    eps_skn = eps_skn.transpose(1, 0, 2).copy()
    occ_skn = occ_skn.transpose(1, 0, 2).copy()
    assert eps_skn.flags.contiguous
    d = dict(nspins=nspins,
             nkpts=nkpts,
             nbands=nbands,
             eigenvalues=eps_skn,
             occupations=occ_skn)
    if fermilevel is not None:
        d.update(efermi=fermilevel)
    return d

def read_static_info_energy(fd, energy_unit):
    def get(name):
        for line in fd:
            if line.strip().startswith(name):
                return float(line.split('=')[-1].strip()) * energy_unit
    return dict(energy=get('Total'), free_energy=get('Free'))


def read_static_info(fd):
    results = {}

    def get_energy_unit(line):  # Convert "title [unit]": ---> unit
        return {'[eV]': eV, '[H]': Hartree}[line.split()[1].rstrip(':')]

    for line in fd:
        if line.strip('*').strip().startswith('Brillouin zone'):
            results.update(read_static_info_kpoints(fd))
        elif line.startswith('Eigenvalues ['):
            unit = get_energy_unit(line)
            results.update(read_static_info_eigenvalues(fd, unit))
        elif line.startswith('Energy ['):
            unit = get_energy_unit(line)
            results.update(read_static_info_energy(fd, unit))
        elif line.startswith('Stress tensor'):
            assert line.split()[-1] == '[H/b^3]'
            stress = read_static_info_stress(fd)
            stress *= Hartree / Bohr**3
            results.update(stress=stress)
        elif line.startswith('Total Magnetic Moment'):
            if 0:
                line = next(fd)
                values = line.split()
                results['magmom'] = float(values[-1])

                line = next(fd)
                assert line.startswith('Local Magnetic Moments')
                line = next(fd)
                assert line.split() == ['Ion', 'mz']
                # Reading  Local Magnetic Moments
                mag_moment = []
                for line in fd:
                    if line == '\n':
                        break  # there is no more thing to search for
                    line = line.replace('\n', ' ')
                    values = line.split()
                    mag_moment.append(float(values[-1]))

                results['magmoms'] = np.array(mag_moment)
        elif line.startswith('Dipole'):
            assert line.split()[-1] == '[Debye]'
            dipole = [float(next(fd).split()[-1]) for i in range(3)]
            results['dipole'] = np.array(dipole) * Debye
        elif line.startswith('Forces'):
            forceunitspec = line.split()[-1]
            forceunit = {'[eV/A]': eV / Angstrom,
                         '[H/b]': Hartree / Bohr}[forceunitspec]
            forces = []
            line = next(fd)
            assert line.strip().startswith('Ion')
            for line in fd:
                if line.strip().startswith('---'):
                    break
                tokens = line.split()[-3:]
                forces.append([float(f) for f in tokens])
            results['forces'] = np.array(forces) * forceunit
        elif line.startswith('Fermi'):
            tokens = line.split()
            unit = {'eV': eV, 'H': Hartree}[tokens[-1]]
            eFermi = float(tokens[-2]) * unit
            results['efermi'] = eFermi

    if 'ibz_k_points' not in results:
        results['ibz_k_points'] = np.zeros((1, 3))
        results['k_point_weights'] = np.ones(1)
    if 0: #'efermi' not in results:
        # Find HOMO level.  Note: This could be a very bad
        # implementation with fractional occupations if the Fermi
        # level was not found otherwise.
        all_energies = results['eigenvalues'].ravel()
        all_occupations = results['occupations'].ravel()
        args = np.argsort(all_energies)
        for arg in args[::-1]:
            if all_occupations[arg] > 0.1:
                break
        eFermi = all_energies[arg]
        results['efermi'] = eFermi

    return results


@reader
def read_octopus(fd, get_kwargs=False):
    kwargs = parse_input_file(fd)

    # input files may contain internal references to other files such
    # as xyz or xsf.  We need to know the directory where the file
    # resides in order to locate those.  If fd is a real file
    # object, it contains the path and we can use it.  Else assume
    # pwd.
    #
    # Maybe this is ugly; maybe it can lead to strange bugs if someone
    # wants a non-standard file-like type.  But it's probably better than
    # failing 'ase gui somedir/inp'
    try:
        fname = fd.name
    except AttributeError:
        directory = None
    else:
        directory = os.path.split(fname)[0]

    atoms, remaining_kwargs = kwargs2atoms(kwargs, directory=directory)
    if get_kwargs:
        return atoms, remaining_kwargs
    else:
        return atoms
