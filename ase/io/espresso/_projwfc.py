from ase.atoms import Atoms
from ase.calculators.espresso import Projwfc
from ._utils import read_fortran_namelist, time_to_float
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np


def read_projwfc_in(fileobj):
    """Parse a projwfc input file, '.pri'

    projwfc inputs are a fortran-namelist format with custom
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

    calc = Projwfc()
    calc.parameters.update(**data['projwfc'])
    atoms = Atoms(calculator=calc)
    atoms.calc.atoms = atoms

    return atoms


def write_projwfc_in(fd, atoms, **kwargs):
    """
    Create an input file for projwfc.

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.

    """

    projwf = ['&projwfc\n']
    for key, value in atoms.calc.parameters.items():
        if value is True:
            projwf.append(f'   {key:16} = .true.\n')
        elif value is False:
            projwf.append(f'   {key:16} = .false.\n')
        elif value is not None:
            if isinstance(value, Path):
                value = str(value)
            # repr format to get quotes around strings
            projwf.append(f'   {key:16} = {value!r:}\n')
    projwf.append('/\n')

    fd.write(''.join(projwf))


def read_projwfc_out(fd):
    """
    Reads projwfc output files

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

    if isinstance(fd, str):
        fd = open(fd, 'r')

    flines = fd.readlines()
    structure = Atoms()

    job_done = False
    walltime = None

    for line in flines:
        if 'JOB DONE' in line:
            job_done = True
        if 'PROJWFC      :' in line:
            time_str = line.split()[-2]
            walltime = time_to_float(time_str)

    calc = Projwfc(atoms=structure)
    calc.results['job done'] = job_done
    calc.results['walltime'] = walltime

    structure.calc = calc

    yield structure

def read_projwfc_xml(fd):
    """
    Reads the atomic_proj.xml file and returns a dictionary of the data contained within it
    """
    tree = ET.parse(fd)
    root = tree.getroot()

    # Reading <HEADER> section
    dct = {k.lower(): (int(v) if '.' not in v else float(v)) for k, v in root[0].attrib.items()}

    # Reading <EIGENSTATES> section
    estates = root[1]
    dct['projectability_data'] = []
    for kpt_element, eigs_element, projs_element in zip(estates[::3], estates[1::3], estates[2::3]):
        kpt = np.array(kpt_element.text.split(), dtype=np.float64)
        eigs = np.array(eigs_element.text.split(), dtype=np.float64)
        projs = []
        for proj_element in projs_element:
            proj = {k.lower(): int(v) for k, v in proj_element.attrib.items()}
            proj['array'] = np.array([float(row.split()[0]) + 1j * float(row.split()[1]) for row in proj_element.text.split('\n') if row.strip()])
            projs.append(proj)
        dct['projectability_data'].append({'kpoint': kpt, 'eigenvalues': eigs, 'projectability': projs})

    return dct
