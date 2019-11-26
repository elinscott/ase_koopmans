import re
import numpy as np

from ase import Atoms
from ase.units import Hartree, Bohr
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)


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
        r'\n[\s]+Geometry \"[\s\S]+\" -> \"[\s\S]*\"[\s]*\n'
        r'^[\s]+[-]+\n\n'
        r'^[\S\s]+\n\n'
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
        r'^[ \t]+[\S]+[ \t]+ENERGY GRADIENTS[ \t]*[\n]+'
        r'^[ \t]+atom[ \t]+coordinates[ \t]+gradient[ \t]*\n'
        r'^(?:[ \t]+x[ \t]+y[ \t]+z){2}[ \t]*\n'
        r'((?:^(?:[ \t]+[\S]+){8}\n)+)\n',
        """\
                         UHF ENERGY GRADIENTS

    atom               coordinates                        gradient
                 x          y          z           x          y          z
   1 C       0.293457  -0.293457   0.293457   -0.000083   0.000083  -0.000083
   2 H       1.125380   1.355351   1.125380    0.000086   0.000089   0.000086
   3 H      -1.355351  -1.125380   1.125380   -0.000089  -0.000086   0.000086
   4 H       1.125380  -1.125380  -1.355351    0.000086  -0.000086  -0.000089
""", re.M)

_nwpw_grad = _define_pattern(
        r'^[ \t]+[=]+[ \t]+Ion Gradients[ \t]+[=]+[ \t]*\n'
        r'^[ \t]+Ion Forces:[ \t]*\n'
        r'((?:^(?:[ \t]+[\S]+){7}\n)+)',
        """\
          =============  Ion Gradients =================
 Ion Forces:
        1 O    (   -0.000012    0.000027   -0.005199 )
        2 H    (    0.000047   -0.013082    0.020790 )
        3 H    (    0.000047    0.012863    0.020786 )
        C.O.M. (   -0.000000   -0.000000   -0.000000 )
          ===============================================
""", re.M)

_paw_grad = _define_pattern(
        r'^[ \t]+[=]+[ \t]+Ion Gradients[ \t]+[=]+[ \t]*\n'
        r'^[ \t]+Ion Positions:[ \t]*\n'
        r'((?:^(?:[ \t]+[\S]+){7}\n)+)'
        r'^[ \t]+Ion Forces:[ \t]*\n'
        r'((?:^(?:[ \t]+[\S]+){7}\n)+)',
        #r'^[\s]+[=]+[\s]+Ion Gradients[\s]+[=]+[\s]*\n'
        #r'^[\s]Ion Positions:[\s]*\n'
        #r'((?:^(?:[\s]+[\S]+){7}\n)+)'
        #r'^[\s]Ion Forces:[\s]*\n'
        #r'((?:^(?:[\s]+[\S]+){7}\n)+)',
        """\
          =============  Ion Gradients =================
 Ion Positions:
        1 O    (   -3.77945   -5.22176   -3.77945 )
        2 H    (   -3.77945   -3.77945    3.77945 )
        3 H    (   -3.77945    3.77945    3.77945 )
 Ion Forces:
        1 O    (   -0.00001   -0.00000    0.00081 )
        2 H    (    0.00005   -0.00026   -0.00322 )
        3 H    (    0.00005    0.00030   -0.00322 )
        C.O.M. (   -0.00000   -0.00000   -0.00000 )
          ===============================================
""", re.M)


_e_gto = dict(mf=_define_pattern(
                    r'^[\s]+Total (?:DFT|SCF) energy =[\s]+([\S]+)[\s]*\n',
                    "         Total SCF energy =    -75.585555997789", re.M),
              mp2=_define_pattern(
                    r'^[\s]+Total MP2 energy[\s]+([\S]+)[\s]*\n',
                    "          Total MP2 energy           -75.708800087578",
                    re.M),
              ccsd=_define_pattern(
                    r'^[\s]+Total CCSD energy:[\s]+([\S]+)[\s]*\n',
                    " Total CCSD energy:            -75.716168566598569",
                    re.M),
              tce=_define_pattern(
                    r'^[\s]+[\S]+[\s]+total energy \/ hartree[\s]+'
                    r'=[\s]+([\S]+)[\s]*\n',
                    " CCD total energy / hartree       "
                    "=       -75.715332545665888", re.M),
              )

_nwpw_energy = _define_pattern(r'^[\s]+Total (?:PSPW|BAND|PAW) energy'
                               r'[\s]+:[\s]+([\S]+)[\s]*\n',
                               " Total PSPW energy     :  -0.1709317826E+02",
                               re.M)

_cell_block = _define_pattern(r'^[ \t]+Lattice Parameters[ \t]*\n'
                              r'^[ \t]+[-]+[ \t]*\n\n[ \t\S]+\n\n'
                              r'((?:^(?:[ \t]+[\S]+){5}\n){3})',
                              """\
      Lattice Parameters
      ------------------

      lattice vectors in angstroms (scale by  1.889725989 to convert to a.u.)

      a1=<   4.000   0.000   0.000 >
      /a2=<   0.000   5.526   0.000 >
      a3=<   0.000   0.000   4.596 >
      a=       4.000 b=      5.526 c=       4.596
      alpha=  90.000 beta=  90.000 gamma=  90.000
      omega=   101.6
""", re.M)

_eval_block = _define_pattern(
        r'^[ \t]+[\S]+ Final (?:Alpha |Beta )?Molecular Orbital Analysis[ \t]*\n'
        r'^[ \t-]+\n\n'
        r'(?:^[ \t]+Vector [ \t\S]+\n(?:^[ \t\S]+\n){3}'
        r'(?:^(?:(?:[ \t]+[\S]+){5}){1,2}[ \t]*\n)+\n)+',
        """\
                       ROHF Final Molecular Orbital Analysis
                       -------------------------------------

 Vector    1  Occ=2.000000D+00  E=-2.043101D+01
              MO Center=  1.1D-20,  1.5D-18,  1.2D-01, r^2= 1.5D-02
   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function  
  ----- ------------  ---------------      ----- ------------  ---------------
     1      0.983233  1 O  s          

 Vector    2  Occ=2.000000D+00  E=-1.324439D+00
              MO Center= -2.1D-18, -8.6D-17, -7.1D-02, r^2= 5.1D-01
   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function  
  ----- ------------  ---------------      ----- ------------  ---------------
     6      0.708998  1 O  s                  1     -0.229426  1 O  s          
     2      0.217752  1 O  s          
""", re.M)

_nwpw_eval_block = _define_pattern(
        r'^[ \t]+(?:virtual )?orbital energies:\n'
        r'(?:^(?:(?:[ \t]+[\S]+){3,4}){1,2}[ \t]*\n)+',
        """\
 orbital energies:
     0.1631366E+00 (   4.439eV)  occ=0.000
     0.1317629E+00 (   3.585eV)  occ=0.000     0.1359902E+00 (   3.701eV)  occ=0.000
     0.1317623E+00 (   3.585eV)  occ=0.000     0.1359899E+00 (   3.701eV)  occ=0.000
    -0.1855813E-01 (  -0.505eV)  occ=0.000    -0.1269164E-01 (  -0.345eV)  occ=0.000
    -0.1804188E+00 (  -4.909eV)  occ=1.000    -0.9891840E-01 (  -2.692eV)  occ=0.000
    -0.3220937E+00 (  -8.765eV)  occ=1.000    -0.3061376E+00 (  -8.330eV)  occ=1.000
    -0.3221004E+00 (  -8.765eV)  occ=1.000    -0.3061414E+00 (  -8.331eV)  occ=1.000
    -0.5795222E+00 ( -15.770eV)  occ=1.000    -0.5516992E+00 ( -15.013eV)  occ=1.000
""", re.M)

_extract_vector = _define_pattern(
        r'^[ \t]+Vector[ \t]+([\S])+[ \t]+Occ=([\S]+)[ \t]+E=([\S]+)[ \t]*\n',
        " Vector    1  Occ=2.000000D+00  E=-2.043101D+01", re.M)

def _get_gto_evals(chunk):
    spin = 1 if re.match(r'[ \t\S]+Beta', chunk) else 0
    data = []
    for vector in _extract_vector.finditer(chunk):
        data.append([float(x.replace('D', 'E')) for x in vector.groups()[1:]])
    data = np.array(data)
    occ = data[:, 0]
    energies = data[:, 1] * Hartree

    return SinglePointKPoint(1., spin, 0, energies, occ)


def _get_gto_kpts(chunk):
    eval_blocks = _eval_block.findall(chunk)
    if not eval_blocks:
        return []
    kpts = []
    kpt = _get_gto_evals(eval_blocks[-1])
    if kpt.s == 1:
        kpts.append(_get_gto_evals(eval_blocks[-2]))
    kpts.append(kpt)
    return kpts


# We support the following properties:
# - Energy (extrapolated, if applicable)
# - Free energy
# - Gradients
# - Dipole moment
# - Quadrupole moment
# - Eigenvalues
# - Occupations

def _parse_geomblock(chunk):
    geomblocks = _geom.findall(chunk)
    if not geomblocks:
        return None
    geomblock = geomblocks[-1].strip().split('\n')
    natoms = len(geomblock)
    symbols = []
    pos = np.zeros((natoms, 3))
    for i, line in enumerate(geomblock):
        line = line.strip().split()
        symbols.append(line[1])
        pos[i] = [float(x) for x in line[3:6]]

    cellblocks = _cell_block.findall(chunk)
    if cellblocks:
        cellblock = cellblocks[-1].strip().split('\n')
        cell = np.zeros((3, 3))
        for i, line in enumerate(cellblock):
            line = line.strip().split()
            cell[i] = [float(x) for x in line[1:4]]
    else:
        cell = None
    return Atoms(symbols, positions=pos, cell=cell)


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

    kpts = _get_gto_kpts(chunk)

    atoms_new = _parse_geomblock(chunk)
    if atoms_new is not None:
        atoms = atoms_new

    if atoms is None:
        return

    calc = SinglePointDFTCalculator(energy=energy, forces=forces, atoms=atoms)
    calc.kpts = kpts
    atoms.set_calculator(calc)
    return atoms


def parse_pw_chunk(chunk):
    atoms = _parse_geomblock(chunk)
    if atoms is None:
        return

    with open('test.txt', 'w') as f:
        f.write(chunk)

    forces = None
    energy = None
    matches = _nwpw_energy.findall(chunk)

    if matches:
        energy = float(matches[-1].replace('D', 'E')) * Hartree

    gradblocks = _nwpw_grad.findall(chunk)
    if gradblocks:
        gradblock = gradblocks[-1].strip().split('\n')
        natoms = len(gradblock)
        symbols = []
        forces = np.zeros((natoms, 3))
        for i, line in enumerate(gradblock):
            line = line.strip().split()
            symbols.append(line[1])
            forces[i] = [float(x) for x in line[3:6]]
        forces *= Hartree / Bohr

    kpts = _get_pw_kpts(chunk)

    calc = SinglePointDFTCalculator(energy=energy, forces=forces, atoms=atoms)
    calc.kpts = kpts
    atoms.set_calculator(calc)
    return atoms



def read_nwchem_out(fobj, index=-1):
    """Splits an NWChem output file into chunks corresponding to
    individual single point calculations."""
    lines = fobj.readlines()

    if index == slice(-1, None, None):
        for line in lines:
            if _gauss_block.match(line):
                return [parse_gto_chunk(''.join(lines))]
            if _pw_block.match(line):
                return [parse_pw_chunk(''.join(lines))]
        else:
            raise ValueError('This does not appear to be a valid NWChem '
                             'output file.')


    # First, find each SCF block
    group = []
    atomslist = []
    header = True
    lastgroup = []
    lastparser = None
    for line in lines:
        group.append(line)
        if _gauss_block.match(line):
            next_parser = parse_gto_chunk
        elif _pw_block.match(line):
            next_parser = parse_pw_chunk
        else:
            continue

        if header:
            header = False
        else:
            atoms = parser(''.join(group))
            if atoms is None and parser is lastparser:
                atoms = parser(''.join(lastgroup + group))
                if atoms is not None:
                    atomslist[-1] = atoms
                    lastgroup += group
            else:
                atomslist.append(atoms)
                lastgroup = group
                lastparser = parser
            group = []
        parser = next_parser
    else:
        if not header:
            atoms = parser(''.join(group))
            if atoms is not None:
                atomslist.append(atoms)

    return atomslist[index]
