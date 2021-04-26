"""Reads Wannier to KC files

Read structures and results from koopmans_ham.x output files. Read
structures from wannier_to_kc.x input files.

"""

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.utils import basestring
from ase.io.espresso._utils import units, Namelist, read_fortran_namelist, generic_construct_namelist
from ase.io.espresso._pw import read_espresso_in, write_espresso_in
from ase.calculators.wann2kc import Wann2KC


KEYS = Namelist((
    ('control', ['prefix', 'outdir', 'kc_iverbosity', 'kc_at_ks', 'homo_only', 'read_unitary_matrix', 'l_vcut']),
    ('wannier', ['seedname', 'check_ks', 'num_wann_occ', 'num_wann_emp', 'have_empty', 'has_disentangle'])))


def write_wann2kc_in(fd, atoms, input_data=None, pseudopotentials=None,
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

    fd.writelines(lines)


def read_wann2kc_in(fileobj):
    data, _ = read_fortran_namelist(fileobj)
    calc = Wann2KC(input_data=data)
    return Atoms(calculator=calc)


def read_wann2kc_out(fileobj):
    """Reads Wannier to KC output files.

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename

    Yields
    ------
    structure : Atoms
        The next structure from the index slice. The Atoms has a
        SinglePointCalculator attached with any results parsed from
        the file.

    """

    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')

    # work with a copy in memory for faster random access
    flines = fileobj.readlines()

    # For the moment, provide an empty atoms object
    structure = Atoms()

    # Extract calculation results
    job_done = False
    for i_line, line in enumerate(flines):
        if 'JOB DONE' in line:
            job_done = True

    # Return an empty calculator object with ths solitary result 'job done'
    calc = SinglePointDFTCalculator(structure)
    calc.results['job_done'] = job_done
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=True, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
