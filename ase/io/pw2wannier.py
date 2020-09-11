"""Reads pw2wannier files.

"""

from ase.atoms import Atoms
from ase.calculators.pw2wannier import PW2Wannier
from ase.io.espresso import Namelist, read_fortran_namelist
import ipdb

def read_pw2wannier_in(fileobj):
    """Parse a pw2wannier input file, '.pw2wan'

    pw2wannier inputs are a fortran-namelist format with custom
    blocks of data. The namelist is parsed as a dict and an atoms object
    is constructed from the included information.

    Parameters
    ----------
    fileobj : file | str
        A file-like object that supports line iteration with the contents
        of the input file, or a filename.

    Returns
    -------
    atoms : Atoms
        Structure defined in the input file.

    Raises
    ------
    KeyError
        Raised for missing keys that are required to process the file
    """
    # TODO: use ase opening mechanisms
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'rU')

    # parse namelist section and extract remaining lines
    data, _ = read_fortran_namelist(fileobj)
    
    if 'inputpp' not in data:
        raise KeyError('Required section &inputpp not found.')

    calc = PW2Wannier()
    calc.parameters['inputpp'] = data['inputpp']
    ipdb.set_trace()
    atoms = Atoms(calculator=calc)
    atoms.calc.atoms = atoms

    return atoms


def write_pw2wannier_in(fd, atoms):
    """
    Create an input file for pw2wannier.

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.

    """

    # Convert to a namelist to make working with parameters much easier
    # Note that the name ``input_data`` is chosen to prevent clash with
    # ``parameters`` in Calculator objects
    if 'inputpp' not in atoms.calc.parameters:
        raise ValueError('No inputpp block found')

    p2w = ['&inputpp\n']
    for key, value in atoms.calc.parameters['inputpp'].items():
        if value is True:
            p2w.append('   {0:16} = .true.\n'.format(key))
        elif value is False:
            p2w.append('   {0:16} = .false.\n'.format(key))
        elif value is not None:
            # repr format to get quotes around strings
            p2w.append('   {0:16} = {1!r:}\n'.format(key, value))
    p2w.append('/\n')

    fd.write(''.join(p2w))
