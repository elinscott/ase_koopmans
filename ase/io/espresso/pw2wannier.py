"""Reads pw2wannier files.

"""

from ase.utils import basestring
from ase.atoms import Atoms
from ase.calculators.espresso import PW2Wannier
from .utils import read_fortran_namelist, time_to_float


def read_pw2wannier_in(fileobj):
    """Parse a pw2wannier input file, 'pw2wan', '.p2wi'

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
    calc.parameters.update(data['inputpp'])
    atoms = Atoms(calculator=calc)
    atoms.calc.atoms = atoms

    return atoms


def write_pw2wannier_in(fd, atoms, **kwargs):
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
    if 'inputpp' not in atoms.calc.parameters:
        raise ValueError('No inputpp block found')

    p2w = ['&inputpp\n']
    for key, value in atoms.calc.parameters.items():
        if value is True:
            p2w.append('   {0:16} = .true.\n'.format(key))
        elif value is False:
            p2w.append('   {0:16} = .false.\n'.format(key))
        elif value is not None:
            # repr format to get quotes around strings
            p2w.append('   {0:16} = {1!r:}\n'.format(key, value))
    p2w.append('/\n')

    fd.write(''.join(p2w))


def read_pw2wannier_out(fd):
    """
    Reads pw2wannier output files

    Parameters
    ----------
    fd : file|str
        A file like object or filename

    Yields
    ------
    structure : atoms
        An Atoms object with an attached SinglePointCalculator containing
        any parsed results
    """

    if isinstance(fd, basestring):
        fd = open(fd, 'rU')

    flines = fd.readlines()

    structure = Atoms()

    job_done = False
    walltime = None

    for line in flines:
        if 'JOB DONE' in line:
            job_done = True
        if 'PW2WANNIER   :' in line:
            time_str = line.split()[-2]
            walltime = time_to_float(time_str)

    calc = PW2Wannier(atoms=structure)
    calc.results['job done'] = job_done
    calc.results['walltime'] = walltime

    structure.calc = calc

    yield structure
