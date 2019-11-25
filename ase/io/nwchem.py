import os
import numpy as np
from ase import Atoms
from ase.units import Hartree, Bohr
import re
from ase.calculators.singlepoint import SinglePointCalculator

_re_symb = r'[A-Z][a-z]?'
_re_int = r'[+-]?[0-9]+'
_re_float = r'[-+]?[0-9]+\.[0-9]+(?:[eEdD][-+][0-9]+)?'

_pattern_test_data = []


def _define_pattern(pattern, example, *args):
    """Accepts a regular expression pattern, as well as an example
    string that the pattern should match. Returns the compiled
    pattern. Additionally, stores the compiled pattern and the
    example string in pattern_test_data for unit testing."""
    regex = re.compile(pattern, *args)
    _pattern_test_data.append((regex, example))
    return regex


_gauss_block = _define_pattern(
        r'^[\s]+NWChem (?:SCF|DFT) Module\n$',
        "                                 NWChem SCF Module\n")

_pw_block = _define_pattern(
        r'^[\s]+\*[\s]+NWPW (?:PSPW|BAND|PAW) Calculation[\s]+\*[\s]*\n$',
        "          *               NWPW PSPW Calculation              *\n")

_geom = _define_pattern(
        r'\n[\s]+Geometry \".+\" -> ".*"[\s]*\n'
        r'^[\s]+[-]+\n'
        r'^[\S\s]+'
        r'^[\s]+No\.[\s]+Tag[\s]+Charge[\s]+X[\s]+Y[\s]+Z\n'
        r'^[-\s]+\n'
        r'((?:^(?:[\s]+[\S]+){6}[\s]*\n)+)',
        """\
                             Geometry "geometry" -> ""
                             -------------------------

 Output coordinates in angstroms (scale by  1.889725989 to convert to a.u.)

  No.       Tag          Charge          X              Y              Z
 ---- ---------------- ---------- -------------- -------------- --------------
    1 C                    6.0000     0.00000000     0.00000000     0.00000000
    2 H                    1.0000     0.62911800     0.62911800     0.62911800
    3 H                    1.0000    -0.62911800    -0.62911800     0.62911800
    4 H                    1.0000     0.62911800    -0.62911800    -0.62911800
""", re.M)

_gto_grad = _define_pattern(
        r'^[\s]+[\S]+[\s]+ENERGY GRADIENTS[\s]*[\n]+'
        r'^[\s]+atom[\s]+coordinates[\s]+gradient[\s]*\n'
        r'^(?:[\s]+x[\s]+y[\s]+z){2}[\s]*\n'
        r'((?:^(?:[\s]+[\S]+){8}\n)+)\n',
        """\
                         UHF ENERGY GRADIENTS

    atom               coordinates                        gradient
                 x          y          z           x          y          z
   1 C       0.293457  -0.293457   0.293457   -0.000083   0.000083  -0.000083
   2 H       1.125380   1.355351   1.125380    0.000086   0.000089   0.000086
   3 H      -1.355351  -1.125380   1.125380   -0.000089  -0.000086   0.000086
   4 H       1.125380  -1.125380  -1.355351    0.000086  -0.000086  -0.000089
""", re.M)

_nwpw_grad = re.compile(r'^[\s]+[=]+[\s]+Ion Gradients[\s]+[=]+[\s]*\n'
                        r'^[\s]+Ion Positions:[\s]*\n'
                        r'^((?:[\s]+{int}\[s]+{symb}\[s]+\(\[s]+'
                        r'(?:{float}){{3}}[\s]+\)[\s]*\n)+)'
                        r'^[\s]+Ion Forces:[\s]*\n'
                        r'^((?:[\s]+{int}\[s]+{symb}\[s]+\(\[s]+'
                        r'(?:{float}){{3}}[\s]+\)[\s]*\n)+)'
                        .format(float=_re_float, int=_re_int,
                                symb=_re_symb), re.M)


_e_gto = dict(mf=re.compile(r'^[\s]+Total (?:DFT|SCF) energy ='
                            r'[\s]+({float})[\s]*\n'.format(float=_re_float),
                            re.M),
              mp2=re.compile(r'^[\s]+Total MP2 energy[\s]+({float})[\s]*\n'
                             r''.format(float=_re_float), re.M),
              ccsd=re.compile(r'^[\s]+Total CCSD energy:[\s]+({float})'
                              r'[\s]*\n'.format(float=_re_float), re.M),
              tce=re.compile(r'^[\s]+[\S]+[\s]+total energy \/ hartree[\s]+'
                             r'=[\s]+({float})[\s]*\n'
                             r''.format(float=_re_float), re.M),
              )

_nwpw_energy = re.compile(r'^[\s]+total[\s]+energy[\s]+:[\s]+({float})'
                          r'[\s]+\([\s]+{float}\/ion\)[\s]*\n'
                          r''.format(float=_re_float), re.M)


# We support the following properties:
# - Energy (extrapolated, if applicable)
# - Free energy
# - Gradients
# - Dipole moment
# - Quadrupole moment
# - Eigenvalues
# - Occupations


def parse_gto_chunk(chunk):
    atoms = None
    forces = None
    energy = None
    for theory in ['tce', 'ccsd', 'mp2', 'mf']:
        matches = _e_gto[theory].findall(chunk)
        if matches:
            energy = float(matches[-1].replace('D', 'E')) * Hartree
            break

    gradblocks = _gto_grad.findall(chunk)
    if gradblocks:
        gradblock = gradblocks[-1].strip().split('\n')
        natoms = len(gradblock)
        symbols = []
        pos = np.zeros((natoms, 3))
        forces = np.zeros((natoms, 3))
        for i, line in enumerate(gradblock):
            line = line.strip().split()
            symbols.append(line[1])
            pos[i] = [float(x) for x in line[2:5]]
            forces[i] = [-float(x) for x in line[5:8]]
        pos *= Bohr
        forces *= Hartree / Bohr
        atoms = Atoms(symbols, positions=pos)

    geomblocks = _geom.findall(chunk)
    if geomblocks:
        geomblock = geomblocks[-1].strip().split('\n')
        natoms = len(geomblock)
        symbols = []
        pos = np.zeros((natoms, 3))
        for i, line in enumerate(geomblock):
            line = line.strip().split()
            symbols.append(line[1])
            pos[i] = [float(x) for x in line[3:6]]
        atoms = Atoms(symbols, positions=pos)

    if atoms is None:
        return

    calc = SinglePointCalculator(energy=energy, forces=forces, atoms=atoms)
    atoms.set_calculator(calc)
    return atoms


def parse_pw_chunk(chunk):
    pass


def read_nwchem_out(fobj, index=-1):
    """Splits an NWChem output file into chunks corresponding to
    individual single point calculations."""

    # First, find each SCF block
    group = []
    chunks = []
    header = True
    for line in fobj:
        if _gauss_block.match(line) or _pw_block.match(line):
            if header:
                header = False
            else:
                chunks.append(''.join(group))
                group = []
        group.append(line)
    else:
        if not header:
            chunks.append(''.join(group))

    if isinstance(index, int):
        return parse_gto_chunk(chunks[int])
    return [parse_gto_chunk(chunk) for chunk in chunks[index]]


_special_kws = ['center', 'autosym', 'autoz', 'theory', 'basis', 'xc', 'task',
                'pseudopotentials', 'set', 'symmetry', 'label', 'geompar',
                'basispar']

_system_type = {1: 'polymer', 2: 'surface', 3: 'crystal'}


def _get_geom(atoms, **params):
    geom_header = ['geometry units angstrom']
    if not params.get('center', False):
        geom_header.append('nocenter')
    if not params.get('autosym', False):
        geom_header.append('noautosym')
    if not params.get('autoz', False):
        geom_header.append('noautoz')
    if 'geompar' in params:
        geom_header.append(params['geompar'])
    geom = [' '.join(geom_header)]

    outpos = atoms.get_positions()
    pbc = atoms.pbc
    if np.any(pbc):
        scpos = atoms.get_scaled_positions()
        for i, pbci in enumerate(pbc):
            if pbci:
                outpos[:, i] = scpos[:, i]
        npbc = pbc.sum()
        cellpars = atoms.get_cell_lengths_and_angles()
        geom.append('  system {} units angstrom'.format(_system_type[npbc]))
        if npbc == 3:
            geom.append('    lattice_vectors')
            for row in atoms.cell:
                geom.append('      {:20.16e} {:20.16e} {:20.16e}'.format(*row))
        else:
            if pbc[0]:
                geom.append('    lat_a {:20.16e}'.format(cellpars[0]))
            if pbc[1]:
                geom.append('    lat_b {:20.16e}'.format(cellpars[1]))
            if pbc[2]:
                geom.append('    lat_c {:20.16e}'.format(cellpars[2]))
            if pbc[1] and pbc[2]:
                geom.append('    alpha {:20.16e}'.format(cellpars[3]))
            if pbc[0] and pbc[2]:
                geom.append('    beta {:20.16e}'.format(cellpars[4]))
            if pbc[1] and pbc[0]:
                geom.append('    gamma {:20.16e}'.format(cellpars[5]))
        geom.append('  end')

    for i, atom in enumerate(atoms):
        geom.append('{:>4} {:20.16e} {:20.16e} {:20.16e}'
                    ''.format(atom.symbol, *outpos[i]))
        symm = params.get('symmetry')
        if symm is not None:
            geom.append('symmetry {}'.format(symm))
    geom.append('end')
    return geom


def _get_basis(**params):
    basis_in = params.get('basis')
    if 'basispar' in params:
        header = 'basis {} noprint'.format(params['basispar'])
    else:
        header = 'basis noprint'
    basis_out = [header]
    if isinstance(basis_in, str):
        basis_out.append('   * library {}'.format(basis_in))
    else:
        for symbol, ibasis in basis_in.items():
            basis_out.append('{:>4} library {}'.format(symbol, ibasis))
    basis_out.append('end')
    return basis_out


_special_keypairs = [('nwpw', 'simulation_cell'),
                     ('nwpw', 'carr-parinello'),
                     ('nwpw', 'brillouin_zone'),
                     ]

def _format_block(key, val, nindent=0):
    prefix = '  ' * nindent
    if val is None:
        return [prefix + key]

    if not isinstance(val, dict):
        return [prefix + '{} {}'.format(key, str(val))]

    out = [prefix + key]
    for subkey, subval in val.items():
        if (key, subkey) in _special_keypairs:
            out += _format_block(subkey, subval, nindent + 1)
        else:
            if isinstance(subval, dict):
                subval = ' '.join([' '.join([a, str(b)]) 
                                             for a, b in subval.items()])
            out.append(prefix + ' '.join([subkey, str(subval)]))
    out.append(prefix + 'end')
    return out


def _get_other(**params):
    out = []
    for kw, block in params.items():
        if kw in _special_kws:
            continue
        out += _format_block(kw, block)
    return out


def _get_set(**params):
    return ['set {} {}'.format(key, val) for key, val in params.items()]


def _get_theory(**params):
    theory = params.get('theory')
    if theory is not None:
        return theory
    nwpw = params.get('nwpw')
    xc = params.get('xc')
    if xc is None:
        if 'tce' in params:
            return 'tce'
        elif 'ccsd' in params:
            return 'ccsd'
        elif 'mp2' in params:
            return 'mp2'
        elif 'scf' in params:
            return 'scf'
        elif nwpw is not None:
            if 'monkhorst-pack' in nwpw or 'brillouin_zone' in nwpw:
                return 'band'
            return 'pspw'
        return 'scf'
    if xc in ['scf', 'dft', 'mp2', 'ccsd', 'tce', 'pspw', 'band', 'paw']:
        return xc
    return 'dft'


_xc_conv = dict(lda='slater pw91lda',
                pbe='xpbe96 cpbe96',
                revpbe='revpbe cpbe96',
                rpbe='rpbe cpbe96',
                pw91='xperdew91 perdew91',
                )


def _update_mult(magmom_tot, **params):
    theory = params['theory']
    if magmom_tot == 0:
        magmom_mult = 1
    else:
        magmom_mult = np.sign(magmom_tot) * (abs(magmom_tot) + 1)
    if 'scf' in params:
        for kw in ['nopen', 'singlet', 'doublet', 'triplet', 'quartet',
                   'quintet', 'sextet', 'septet', 'octet']:
            if kw in params['scf']:
                break
        else:
            params['scf']['nopen'] = magmom_tot
    elif theory in ['scf', 'mp2', 'ccsd', 'tce']:
        params['scf'] = dict(nopen=magmom_tot)

    if 'dft' in params:
        if 'mult' not in params['dft']:
            params['dft']['mult'] = magmom_mult
    elif theory == 'dft':
        params['dft'] = dict(mult=magmom_mult)

    if 'nwpw' in params:
        if 'mult' not in params['nwpw']:
            params['nwpw']['mult'] = magmom_mult
    elif theory in ['pspw', 'band', 'paw']:
        params['nwpw'] = dict(mult=magmom_mult)

    return params


def write_nwchem_in(fd, atoms, properties=None, **params):
    params = params.copy()
    perm = os.path.abspath(params.get('perm', params['label']))
    scratch = os.path.abspath(params.get('scratch', params['label']))
    os.makedirs(perm, exist_ok=True)
    os.makedirs(scratch, exist_ok=True)

    assert properties is not None
    if properties is None:
        properties = ['energy']

    task = params.get('task')
    if task is None:
        if 'forces' in properties:
            task = 'gradient'
        else:
            task = 'energy'

    theory = _get_theory(**params)
    params['theory'] = theory
    xc = params.get('xc')
    if 'xc' in params:
        xc = _xc_conv.get(params['xc'], params['xc'])
        if theory == 'dft':
            if 'dft' not in params:
                params['dft'] = dict()
            params['dft']['xc'] = xc
        elif theory in ['pspw', 'band', 'paw']:
            if 'nwpw' not in params:
                params['nwpw'] = dict()
            params['nwpw']['xc'] = xc

    magmom_tot = int(atoms.get_initial_magnetic_moments().sum())
    params = _update_mult(magmom_tot, **params)

    out = ['title "{}"'.format(params['label']),
           'permanent_dir {}'.format(perm),
           'scratch_dir {}'.format(scratch),
           'start {}'.format(params['label']),
           '\n'.join(_get_geom(atoms, **params)),
           '\n'.join(_get_basis(**params)),
           '\n'.join(_get_other(**params)),
           '\n'.join(_get_set(**params.get('set', dict()))),
           'task {} {}'.format(theory, task)]

    fd.write('\n\n'.join(out))
