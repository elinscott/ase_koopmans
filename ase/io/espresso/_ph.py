"""Reads pw2wannier/wann2kcp files.

"""

from ase import Atom
from pathlib import Path
from ase.utils import basestring
from ase.atoms import Atoms
from ._utils import read_fortran_namelist, time_to_float
from ase.calculators.espresso import EspressoPh


def read_ph_in(fileobj):
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

    calc = EspressoPh()
    calc.parameters.update(**data['inputph'])
    atoms = Atoms(calculator=calc)
    atoms.calc.atoms = atoms

    return atoms


def write_ph_in(fd, atoms, **kwargs):
    """
    Create an input file for ph.x

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.

    """

    ph = ['&inputph\n']

    masses = {}
    for i, element in enumerate(atoms.calc.parameters.pseudopotentials.keys()):
        masses[f'amass({i+1})'] = Atom(element).mass

    all_parameters = dict(**atoms.calc.parameters, **masses)
    all_parameters.pop('pseudopotentials', None)

    for key, value in all_parameters.items():
        if value is True:
            ph.append('   {0:16} = .true.\n'.format(key))
        elif value is False:
            ph.append('   {0:16} = .false.\n'.format(key))
        elif value is not None:
            if isinstance(value, Path):
                value = str(value)
            # repr format to get quotes around strings
            ph.append('   {0:16} = {1!r:}\n'.format(key, value))
    ph.append('/\n')
    ph.append('0.0 0.0 0.0')

    fd.write(''.join(ph))


def read_ph_out(fd, *args, **kwargs):
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
        if line.strip().startswith('PHONON'):
            time_str = line.split()[-2]
            walltime = time_to_float(time_str)

    calc = EspressoPh(atoms=structure)
    calc.results['job done'] = job_done
    calc.results['walltime'] = walltime

    structure.calc = calc

    yield structure
