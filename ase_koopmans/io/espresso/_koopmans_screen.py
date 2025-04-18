"""Reads Koopmans Screen files

Read structures and results from kcw.x (screen mode) output files.
Read structures from kcw.x (screen mode) input files.
"""

import copy
from ase_koopmans import Atoms
from ase_koopmans.calculators.singlepoint import SinglePointDFTCalculator
from ase_koopmans.utils import base_koopmansstring
from ._utils import read_fortran_namelist, generic_construct_namelist, time_to_float, units, dict_to_input_lines
from ._wann2kc import KEYS as W2KCW_KEYS
from ase_koopmans.calculators.espresso import KoopmansScreen

KEYS = copy.deepcopy(W2KCW_KEYS)
KEYS['SCREEN'] = ['tr2', 'nmix', 'niter', 'eps_inf', 'i_orb', 'check_spread']


def write_koopmans_screen_in(fd, atoms, input_data=None, **kwargs):

    if 'input_data' in atoms.calc.parameters and input_data is None:
        input_data = atoms.calc.parameters['input_data']

    input_parameters = construct_namelist(input_data, **kwargs)

    assert input_parameters['CONTROL']['calculation'] == 'screen'

    lines = []
    for section in input_parameters:
        assert section in KEYS.keys()

        if section == 'WANNIER' and not input_parameters['CONTROL'].get('kcw_at_ks', True):
            # Do not write the WANNIER section if kcw_at_ks is true
            continue

        lines.append('&{0}\n'.format(section.upper()))
        lines += dict_to_input_lines(input_parameters[section])
        lines.append('/\n')  # terminate section

    fd.writelines(lines)


def read_koopmans_screen_in(fileobj):
    data, _ = read_fortran_namelist(fileobj)
    calc = KoopmansScreen(**{k: v for block in data.values() for k, v in block.items()})
    return Atoms(calculator=calc)


def read_koopmans_screen_out(fileobj):
    """Reads Koopmans Screen output files.

    Will probably raise errors for broken or incomplete files.

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
    job_done = False
    walltime = None
    alphas = [[]]
    orbital_data = {'self-Hartree': []}
    for i_line, line in enumerate(flines):
        if 'relaxed' in line:
            splitline = line.split()
            alphas[-1].append(float(splitline[-5]))
            orbital_data['self-Hartree'].append(float(splitline[-1]) * units.Ry)

        if 'JOB DONE' in line:
            job_done = True

        if 'KCW          :' in line:
            time_str = line.split('CPU')[-1].split('WALL')[0].strip()
            walltime = time_to_float(time_str)

    # Put everything together
    calc = SinglePointDFTCalculator(structure)
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    calc.results['alphas'] = alphas
    calc.results['orbital_data'] = orbital_data
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=False, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
