"""Reads pw2wannier/wann2kcp files.

"""

from pathlib import Path
from ase.utils import basestring
from ase.atoms import Atoms
from ._utils import read_fortran_namelist, time_to_float


def read_x2y_in(fileobj, calc_class):
    """Parse a pw2wannier/wann2kcp input file

    inputs are a fortran-namelist format with custom blocks of data.
    The namelist is parsed as a dict and an atoms object is constructed
    from the included information.

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

    calc = calc_class()
    calc.parameters.update(**data['inputpp'])
    atoms = Atoms(calculator=calc)
    atoms.calc.atoms = atoms

    return atoms


def write_x2y_in(fd, atoms, **kwargs):
    """
    Create an input file for pw2wannier/wann2kcp.

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.

    """

    x2y = ['&inputpp\n']
    for key, value in atoms.calc.parameters.items():
        if value is True:
            x2y.append(f'   {key:16} = .true.\n')
        elif value is False:
            x2y.append(f'   {key:16} = .false.\n')
        elif value is not None:
            if isinstance(value, Path):
                value = str(value)
            # repr format to get quotes around strings
            x2y.append(f'   {key:16} = {value!r:}\n')
    x2y.append('/\n')

    fd.write(''.join(x2y))


def read_x2y_out(fd, calc_class):
    """
    Reads pw2wannier/wann2kcp output files

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
        if line.startswith(calc_class.__name__.upper()):
            time_str = line.split()[-2]
            walltime = time_to_float(time_str)

    calc = calc_class(atoms=structure)
    calc.results['job done'] = job_done
    calc.results['walltime'] = walltime

    structure.calc = calc

    yield structure
