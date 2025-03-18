"""Reads pw2wannier/wann2kcp files.

"""

from pathlib import Path

from ase_koopmans import Atom
from ase_koopmans.atoms import Atoms
from ase_koopmans.calculators.espresso import EspressoPh
from ase_koopmans.utils import base_koopmansstring

from ._utils import dict_to_input_lines, read_fortran_namelist, time_to_float


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
    # TODO: use ase_koopmans opening mechanisms
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

    all_parameters = dict(**atoms.calc.parameters)
    all_parameters.pop('pseudopotentials', None)

    ph += dict_to_input_lines(all_parameters)
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

    if isinstance(fd, base_koopmansstring):
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
