"""Reads Wannier to KCW files

Read structures and results from kcw.x (wann2kcw mode) output files.
Read structures from kcw.x (wann2kcw mode) input files.

"""

from ase_koopmans.atoms import Atoms
from ase_koopmans.calculators.espresso import Wann2KC
from ase_koopmans.calculators.singlepoint import SinglePointDFTCalculator
from ase_koopmans.utils import base_koopmansstring

from ._utils import (Namelist, dict_to_input_lines, generic_construct_namelist,
                     read_fortran_namelist, time_to_float)

KEYS = Namelist((
    ('control', ['prefix', 'outdir', 'kcw_iverbosity', 'kcw_at_ks', 'calculation', 'lrpa',
                 'mp1', 'mp2', 'mp3', 'homo_only', 'read_unitary_matrix', 'l_vcut', 'assume_isolated',
                 'spin_component']),
    ('wannier', ['seedname', 'check_ks', 'num_wann_occ', 'num_wann_emp', 'have_empty', 'has_disentangle', 'l_unique_manifold'])))


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
        lines += dict_to_input_lines(input_parameters[section])
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

    if isinstance(fileobj, base_koopmansstring):
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
        if 'KCW          :' in line:
            time_str = line.split('CPU')[-1].split('WALL')[0].strip()
            walltime = time_to_float(time_str)

    # Return an empty calculator object with ths solitary result 'job done'
    calc = SinglePointDFTCalculator(structure)
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=True, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
