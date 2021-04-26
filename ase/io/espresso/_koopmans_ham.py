"""Reads Koopmans Ham files

Read structures and results from koopmans_ham.x output files. Read
structures from koopmans_ham.x input files.
"""

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.utils import basestring
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure
import ase.io.espresso._utils as utils
from ase.io.wann2kc import KEYS as W2KKEYS

from ase.calculators.koopmans_ham import KoopmansHam

KEYS = copy.deepcopy(W2KKEYS)
KEYS['HAM'] = ['do_bands', 'use_ws_distance', 'write_hr', 'l_alpha_corr', 'lrpa', 'mp1', 'mp2', 'mp3']


def write_koopmans_ham_in(fd, atoms, input_data=None, pseudopotentials=None,
                          kspacing=None, kpts=None, koffset=(0, 0, 0), **kwargs):

    if 'input_data' in atoms.calc.parameters and input_data is None:
        input_data = atoms.calc.parameters['input_data']

    input_parameters = construct_namelist(input_data, **kwargs)
    lines = []
    for section in input_parameters:
        assert section in KEYS.keys()
        lines.append('&{0}\n'.format(section.upper()))
        for key, value in input_parameters[section].items():
            if value is True:
                lines.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                lines.append('   {0:16} = .false.\n'.format(key))
            elif value is not None:
                # repr format to get quotes around strings
                lines.append('   {0:16} = {1!r:}\n'.format(key, value))
        lines.append('/\n')  # terminate section

    # kpoints block
    lines += construct_kpoints_card(atoms, kpts, kspacing, koffset)

    fd.writelines(lines)


def read_koopmans_ham_in(fileobj):
    data, _ = read_fortran_namelist(fileobj)
    calc = KoopmansHam(input_data=data)
    return Atoms(calculator=calc)


def read_koopmans_ham_out(fileobj):
    """Reads Koopmans Ham output files.

    Will probably raise errors for broken or incomplete files.

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename

    Yields
    ------
    structure : Atoms
        The next structure. The Atoms has a SinglePointCalculator attached with any results parsed from the file.

    """

    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')

    # work with a copy in memory for faster random access
    flines = fileobj.readlines()

    # For the moment, provide an empty atoms object
    structure = Atoms()

    # Extract calculation results
    job_done = False
    kpts = []
    energies = []
    for i_line, line in enumerate(flines):
        if 'KC interpolated eigenvalues at k=' in line:
            kpts.append([float(x) for x in line.split()[-3:]])
            energies.append([float(x) for x in flines[i_line + 2].split()])

        if 'JOB DONE' in line:
            job_done = True

    # Put everything together
    calc = SinglePointDFTCalculator(structure)
    calc.results['job_done'] = job_done
    calc.results['energies'] = energies
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=False, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
