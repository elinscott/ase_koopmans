import re
from ase.utils import basestring
from ase.units import Hartree, Bohr
import numpy as np
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
        # XXX This processing should probably happen higher up
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


def read_abinit_out(fd):
    text = fd.read().lower()
    results = {}

    for line in iter(text.split('\n')):
        if line.rfind('error') > -1 or line.rfind('was not enough scf cycles to converge') > -1:
            raise ReadError(line)
        if line.rfind('natom  ') > -1:
            natoms = int(line.split()[-1])

    lines = iter(text.split('\n'))
    # Stress:
    # Printed in the output in the following format [Hartree/Bohr^3]:
    # sigma(1 1)=  4.02063464E-04  sigma(3 2)=  0.00000000E+00
    # sigma(2 2)=  4.02063464E-04  sigma(3 1)=  0.00000000E+00
    # sigma(3 3)=  4.02063464E-04  sigma(2 1)=  0.00000000E+00
    for line in lines:
        if line.rfind(
            'cartesian components of stress tensor (hartree/bohr^3)') > -1:
            stress = np.empty(6)
            for i in range(3):
                entries = next(lines).split()
                stress[i] = float(entries[2])
                stress[i + 3] = float(entries[5])
            results['stress'] = stress * Hartree / Bohr**3
            break
    else:
        raise RuntimeError

    # Energy [Hartree]:
    # Warning: Etotal could mean both electronic energy and free energy!
    etotal = None
    efree = None
    if 'PAW method is used'.lower() in text:  # read DC energy according to M. Torrent
        for line in iter(text.split('\n')):
            if line.rfind('>>>>> internal e=') > -1:
                etotal = float(line.split('=')[-1])*Hartree  # second occurrence!
        for line in iter(text.split('\n')):
            if line.rfind('>>>> etotal (dc)=') > -1:
                efree = float(line.split('=')[-1])*Hartree
    else:
        for line in iter(text.split('\n')):
            if line.rfind('>>>>> internal e=') > -1:
                etotal = float(line.split('=')[-1])*Hartree  # first occurrence!
                break
        for line in iter(text.split('\n')):
            if line.rfind('>>>>>>>>> etotal=') > -1:
                efree = float(line.split('=')[-1])*Hartree
    if efree is None:
        raise RuntimeError('Total energy not found')
    if etotal is None:
        etotal = efree

    # Energy extrapolated to zero Kelvin:
    results['energy'] = (etotal + efree) / 2
    results['free_energy'] = efree

    # Forces:
    for line in lines:
        if line.rfind('cartesian forces (ev/angstrom) at end:') > -1:
            forces = []
            for i in range(natoms):
                forces.append(np.array(
                        [float(f) for f in next(lines).split()[1:]]))
            results['forces'] = np.array(forces)
            break
    else:
        raise RuntimeError

    results['nbands'] = read_number_of_bands(lines)
    results['niter'] = read_number_of_iterations(lines)
    results['magmom'] = read_magnetic_moment(lines)
    return results

def read_number_of_iterations(lines):
    niter = None
    for line in lines:
        if ' At SCF step' in line:
            niter = int(line.split()[3].rstrip(','))
    return niter

def read_number_of_bands(lines):
    nband = None
    for line in lines:
        if '     nband' in line: # nband, or nband1, nband*
            nband = int(line.split()[-1].strip())
    return nband

def read_magnetic_moment(lines):
    magmom = 0.0
    # only for spinpolarized system Magnetisation is printed
    for line in lines:
        if 'Magnetisation' in line:
            magmom = float(line.split('=')[-1].strip())
    return magmom


def read_eig(fd):
    line = next(fd)
    results = {}
    m = re.match(r'\s*Fermi \(or HOMO\) energy \(hartree\)\s*=\s*(\S+)',
                 line)
    assert m is not None
    results['fermilevel'] = float(m.group(1)) * Hartree
    # XXX TODO spin
    kpoint_weights = []
    kpoint_coords = []

    eig_skn = []
    for ispin in range(2):
        print('ispin', ispin)
        # (We don't know if we have two spins until we see next line)
        line = next(fd)
        m = re.match(r'\s*Eigenvalues \(hartree\) for nkpt\s*='
                     r'\s*(\S+)\s*k\s*points', line)
        nspins = 2 if 'SPIN' in line else 1
        assert m is not None
        nkpts = int(m.group(1))

        headerpattern = (r'\s*kpt#\s*\S+\s*'
                         r'nband=\s*(\d+),\s*'
                         r'wtk=([^,]+),\s*'
                         r'kpt=\s*(\S)+\s*(\S+)\s*(\S+)')


        eig_kn = []
        for ikpt in range(nkpts):
            header = next(fd)
            m = re.match(headerpattern, header)
            assert m is not None, header
            nbands = int(m.group(1))
            weight = float(m.group(2))
            kvector = np.array(m.group(3, 4, 5)).astype(float)
            if ispin == 0:
                kpoint_coords.append(kvector)
                kpoint_weights.append(float(weight))

            eig_n = []
            while len(eig_n) < nbands:
                line = next(fd)
                tokens = line.split()
                values = np.array(tokens).astype(float) * Hartree
                eig_n.extend(values)
            assert len(eig_n) == nbands
            eig_kn.append(eig_n)
            assert nbands == len(eig_kn[0])
        eig_skn.append(eig_kn)

        if nspins == 1:
            break

    eig_skn = np.array(eig_skn)
    assert eig_skn.shape == (nspins, nkpts, nbands), (eig_skn.shape, (nspins, nkpts, nbands))
    results['ibz_kpoints'] = np.array(kpoint_coords)
    results['kpoint_weights'] = np.array(kpoint_weights)
    results['eigenvalues'] = eig_skn
    return results


def read_abinit_log(fd):
    results = {}
    width = None
    nelect = None
    lines = fd.readlines()
    for line in lines:
        if 'tsmear' in line:
            results['width'] = float(line.split()[1].strip()) * Hartree
        if 'with nelect' in line:
            results['nelect'] = float(line.split('=')[1].strip())
    return results
