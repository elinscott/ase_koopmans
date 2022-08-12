"""Reads Wannier to KCW files

Read structures and results from kcw.x (wann2kcw mode) output files.
Read structures from kcw.x (wann2kcw mode) input files.

"""

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.utils import basestring
from ._utils import Namelist, read_fortran_namelist, generic_construct_namelist, time_to_float
from ase.calculators.espresso import Wann2KC


KEYS = Namelist((
    ('control', ['prefix', 'outdir', 'kcw_iverbosity', 'kcw_at_ks', 'calculation', 'lrpa',
                 'mp1', 'mp2', 'mp3', 'homo_only', 'read_unitary_matrix', 'l_vcut', 'assume_isolated', 'spin_component']),
    ('wannier', ['seedname', 'check_ks', 'num_wann_occ', 'num_wann_emp', 'have_empty', 'has_disentangle'])))


def write_wann2kc_in(fd, atoms, input_data=None, pseudopotentials=None,
                     kspacing=None, kpts=None, koffset=(0, 0, 0), **kwargs):

    if 'input_data' in atoms.calc.parameters and input_data is None:
        input_data = atoms.calc.parameters['input_data']

    input_parameters = construct_namelist(input_data, **kwargs)

    assert input_parameters['CONTROL']['calculation'] == 'wann2kcw'

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
    calc = Wann2KC(**{k: v for block in data.values() for k, v in block.items()})
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
    walltime = None
    job_done = False
    for line in flines:
        if 'JOB DONE' in line:
            job_done = True
        if 'KC_WANN      :' in line:
            time_str = line.split()[-2]
            walltime = time_to_float(time_str)

    # Return an empty calculator object with ths solitary result 'job done'
    calc = SinglePointDFTCalculator(structure)
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=True, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
