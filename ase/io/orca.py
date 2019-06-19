from ase.utils import StringIO
from ase.io import read

# Made from NWChem interface

def read_orca_output(filename):
    """Method to read geometry from a ORCA output."""

    f = filename
    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].find('XYZ format geometry') >= 0:
            natoms = int(lines[i + 2].split()[0])
            string = ''
            for j in range(2, natoms + 4):
                xyzstring = lines[i + j]
                symbol = xyzstring.split()[0].strip()
                # replace bq ghost with X: MDTMP can we do better?
                if symbol.startswith('bq'):
                    xyzstring = xyzstring.replace(symbol, 'X')
                string += xyzstring
            atoms = read(StringIO(string), format='xyz')
            i += natoms + 4
        else:
            i += 1

    if isinstance(filename, str):
        f.close()

    return atoms


def read_orca(filename):
    """Method to read geometry from an ORCA input file."""
    f = filename
    if isinstance(filename, str):
        f = open(filename)
    lines = f.readlines()

    # Find geometry region of input file.
    stopline = 0
    for index, line in enumerate(lines):
        if line.startswith('geometry'):
            startline = index + 1
            stopline = -1
        elif (line.startswith('end') and stopline == -1):
            stopline = index
    # Format and send to read_xyz.
    xyz_text = '%i\n' % (stopline - startline)
    xyz_text += ' geometry\n'
    for line in lines[startline:stopline]:
        xyz_text += line
    atoms = read(StringIO(xyz_text), format='xyz')
    atoms.set_cell((0., 0., 0.))  # no unit cell defined

    if type(filename) == str:
        f.close()

    return atoms


def write_orca(filename, atoms, charge, mult):
    """Function to write ORCA input file coordinates in simple format
    """

    if isinstance(filename, str):
        f = open(filename, 'w')
    else:  # Assume it's a 'file-like object'
        f = filename
    f.write('*xyz');f.write(" %d" % charge);f.write(" %d \n" % mult)
    for atom in atoms:
        if atom.tag == 71:  # 71 is ascii G (Ghost)
            symbol = atom.symbol + ' : '
        else:
            symbol = atom.symbol + '   '
        f.write(symbol +
                str(atom.position[0]) + ' ' +
                str(atom.position[1]) + ' ' +
                str(atom.position[2]) + '\n')
    f.write('*\n')
