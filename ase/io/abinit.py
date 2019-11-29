from ase.utils import basestring
"""
This module contains functionality for reading an ASE
Atoms object in ABINIT input format.

"""

def read_abinit(filename='abinit.in'):
    """Import ABINIT input file.

    Reads cell, atom positions, etc. from abinit input file
    """

    from ase import Atoms, units

    if isinstance(filename, basestring):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    lines = []
    for line in f.readlines():
        meat = line.split('#', 1)[0]
        lines.append(meat)
    tokens = ' '.join(lines).lower().split()

    if isinstance(filename, basestring):
        f.close()

    # note that the file can not be scanned sequentially

    index = tokens.index("acell")
    unit = 1.0
    if(tokens[index + 4].lower()[:3] != 'ang'):
        unit = units.Bohr
    acell = [unit * float(tokens[index + 1]),
             unit * float(tokens[index + 2]),
             unit * float(tokens[index + 3])]

    index = tokens.index("natom")
    natom = int(tokens[index+1])

    index = tokens.index("ntypat")
    ntypat = int(tokens[index+1])

    index = tokens.index("typat")
    typat = []
    for i in range(natom):
        t = tokens[index+1+i]
        if '*' in t:  # e.g. typat 4*1 3*2 ...
            typat.extend([int(t) for t in ((t.split('*')[1] + ' ') * int(t.split('*')[0])).split()])
        else:
            typat.append(int(t))
        if len(typat) == natom: break

    index = tokens.index("znucl")
    znucl = []
    for i in range(ntypat):
        znucl.append(int(tokens[index+1+i]))

    index = tokens.index("rprim")
    rprim = []
    for i in range(3):
        rprim.append([acell[i]*float(tokens[index+3*i+1]),
                      acell[i]*float(tokens[index+3*i+2]),
                      acell[i]*float(tokens[index+3*i+3])])

    # create a list with the atomic numbers
    numbers = []
    for i in range(natom):
        ii = typat[i] - 1
        numbers.append(znucl[ii])

    # now the positions of the atoms
    if "xred" in tokens:
        index = tokens.index("xred")
        xred = []
        for i in range(natom):
            xred.append([float(tokens[index+3*i+1]),
                         float(tokens[index+3*i+2]),
                         float(tokens[index+3*i+3])])
        atoms = Atoms(cell=rprim, scaled_positions=xred, numbers=numbers,
                      pbc=True)
    else:
        if "xcart" in tokens:
            index = tokens.index("xcart")
            unit = units.Bohr
        elif "xangst" in tokens:
            unit = 1.0
            index = tokens.index("xangst")
        else:
            raise IOError(
                "No xred, xcart, or xangs keyword in abinit input file")

        xangs = []
        for i in range(natom):
            xangs.append([unit*float(tokens[index+3*i+1]),
                          unit*float(tokens[index+3*i+2]),
                          unit*float(tokens[index+3*i+3])])
        atoms = Atoms(cell=rprim, positions=xangs, numbers=numbers, pbc=True)
    
    try:
        ii = tokens.index('nsppol')
    except ValueError:
        nsppol = None
    else:
        nsppol = int(tokens[ii + 1])

    if nsppol == 2:
        index = tokens.index('spinat')
        magmoms = [float(tokens[index + 3 * i + 3]) for i in range(natom)]
        atoms.set_initial_magnetic_moments(magmoms)

    return atoms


keys_with_units = {
    'toldfe': 'eV',
    'tsmear': 'eV',
    'paoenergyshift': 'eV',
    'zmunitslength': 'Bohr',
    'zmunitsangle': 'rad',
    'zmforcetollength': 'eV/Ang',
    'zmforcetolangle': 'eV/rad',
    'zmmaxdispllength': 'Ang',
    'zmmaxdisplangle': 'rad',
    'ecut': 'eV',
    'pawecutdg': 'eV',
    'dmenergytolerance': 'eV',
    'electronictemperature': 'eV',
    'oneta': 'eV',
    'onetaalpha': 'eV',
    'onetabeta': 'eV',
    'onrclwf': 'Ang',
    'onchemicalpotentialrc': 'Ang',
    'onchemicalpotentialtemperature': 'eV',
    'mdmaxcgdispl': 'Ang',
    'mdmaxforcetol': 'eV/Ang',
    'mdmaxstresstol': 'eV/Ang**3',
    'mdlengthtimestep': 'fs',
    'mdinitialtemperature': 'eV',
    'mdtargettemperature': 'eV',
    'mdtargetpressure': 'eV/Ang**3',
    'mdnosemass': 'eV*fs**2',
    'mdparrinellorahmanmass': 'eV*fs**2',
    'mdtaurelax': 'fs',
    'mdbulkmodulus': 'eV/Ang**3',
    'mdfcdispl': 'Ang',
    'warningminimumatomicdistance': 'Ang',
    'rcspatial': 'Ang',
    'kgridcutoff': 'Ang',
    'latticeconstant': 'Ang'}


import numpy as np
from ase.calculators.calculator import kpts2mp
def write_abinit_in(fd, atoms, param, species):
    inp = {}
    inp.update(param)
    for key in ['xc', 'smearing', 'kpts', 'pps', 'raw']:
        del inp[key]

    smearing = param.get('smearing')
    if 'tsmear' in param or 'occopt' in param:
        assert smearing is None

    if smearing is not None:
        inp['occopt'] = {'fermi-dirac': 3,
                         'gaussian': 7}[smearing[0].lower()]
        inp['tsmear'] = smearing[1]

    inp['natom'] = len(atoms)

    if 'nbands' in param:
        inp['nband'] = param.nbands
        del inp['nbands']

    # ixc is set from paw/xml file. Ignore 'xc' setting then.
    if param.get('pps') not in ['pawxml']:
        if 'ixc' not in param:
            inp['ixc'] = {'LDA': 7,
                          'PBE': 11,
                          'revPBE': 14,
                          'RPBE': 15,
                          'WC': 23}[param.xc]

    magmoms = atoms.get_initial_magnetic_moments()
    if magmoms.any():
        inp['nsppol'] = 2
        fd.write('spinat\n')
        for n, M in enumerate(magmoms):
            fd.write('%.14f %.14f %.14f\n' % (0, 0, M))
    else:
        inp['nsppol'] = 1

    for key in sorted(inp.keys()):
        value = inp[key]
        unit = keys_with_units.get(key)
        if unit is None:
            fd.write('%s %s\n' % (key, value))
        else:
            if 'fs**2' in unit:
                value /= fs**2
            elif 'fs' in unit:
                value /= fs
            fd.write('%s %e %s\n' % (key, value, unit))

    if param.raw is not None:
        for line in param.raw:
            if isinstance(line, tuple):
                fd.write(' '.join(['%s' % x for x in line]) + '\n')
            else:
                fd.write('%s\n' % line)

    fd.write('#Definition of the unit cell\n')
    fd.write('acell\n')
    fd.write('%.14f %.14f %.14f Angstrom\n' % (1.0, 1.0, 1.0))
    fd.write('rprim\n')
    if atoms.number_of_lattice_vectors != 3:
        raise RuntimeError('Abinit requires a 3D cell, but cell is {}'
                           .format(atoms.cell))
    for v in atoms.cell:
        fd.write('%.14f %.14f %.14f\n' %  tuple(v))

    fd.write('chkprim 0 # Allow non-primitive cells\n')

    fd.write('#Definition of the atom types\n')
    fd.write('ntypat %d\n' % (len(species)))
    fd.write('znucl')
    for n, Z in enumerate(species):
        fd.write(' %d' % (Z))
    fd.write('\n')
    fd.write('#Enumerate different atomic species\n')
    fd.write('typat')
    fd.write('\n')
    types = []
    for Z in atoms.numbers:
        for n, Zs in enumerate(species):
            if Z == Zs:
                types.append(n + 1)
    n_entries_int = 20  # integer entries per line
    for n, type in enumerate(types):
        fd.write(' %d' % (type))
        if n > 1 and ((n % n_entries_int) == 1):
            fd.write('\n')
    fd.write('\n')

    fd.write('#Definition of the atoms\n')
    fd.write('xangst\n')
    for pos in atoms.positions:
        fd.write('%.14f %.14f %.14f\n' %  tuple(pos))

    if 'kptopt' not in param:
        mp = kpts2mp(atoms, param.kpts)
        fd.write('kptopt 1\n')
        fd.write('ngkpt %d %d %d\n' % tuple(mp))
        fd.write('nshiftk 1\n')
        fd.write('shiftk\n')
        fd.write('%.1f %.1f %.1f\n' % tuple((np.array(mp) + 1) % 2 * 0.5))

    fd.write('chkexit 1 # abinit.exit file in the running directory terminates after the current SCF\n')


def write_abinit(filename, atoms, cartesian=False, long_format=True):
    """Method to write abinit input files."""

    import numpy as np
    from ase import data

    if isinstance(filename, basestring):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename

    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                            "one image to input")
        else:
            atoms = atoms[0]

    # Write atom positions in scaled or cartesian coordinates
    if cartesian:
        coord = atoms.get_positions()
    else:
        coord = atoms.get_scaled_positions()

    # let us order the atoms according to chemical symbol
    ind = np.argsort(atoms.get_chemical_symbols())
    symbols = np.array(atoms.get_chemical_symbols())[ind]
    coord = coord[ind]

    # and now we count how many atoms of which type we have
    sc = []
    psym = symbols[0]
    count = 0
    for sym in symbols:
        if sym != psym:
            sc.append((psym, count))
            psym = sym
            count = 1
        else:
            count += 1
    sc.append((psym, count))

    f.write('\n# Definition of the atom types\n')
    f.write("ntypat  " + str(len(sc)) + "\n")
    f.write("znucl  ")
    for specie in sc:
        f.write(str(data.atomic_numbers[specie[0]]) + " ")
    f.write('\n')

    f.write('\n# Definition of the atoms\n')
    f.write('natom  ' + str(len(symbols)) + '\n')
    f.write('typat  ')
    typat = 1
    for specie in sc:
        for natom in range(specie[1]):
            f.write(str(typat) + ' ')
        typat = typat + 1
    f.write('\n')

    f.write('\n# Definition of the unit cell\n')
    f.write('acell\n')
    f.write('%.14f %.14f %.14f Angstrom\n' %  (1.0, 1.0, 1.0))
    f.write('\n')
    f.write('rprim\n')
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'

    for vec in atoms.get_cell():
        f.write(' ')
        for el in vec:
            f.write(latt_form % el)
        f.write('\n')
    f.write('\n')

    # Write atom positions in scaled or cartesian coordinates
    if cartesian:
        f.write('xangst\n')
    else:
        f.write('xred\n')

    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'

    for iatom, atom in enumerate(coord):
        f.write(' ')
        for dcoord in atom:
            f.write(cform % dcoord)
        f.write('\n')

    if isinstance(filename, basestring):
        f.close()
