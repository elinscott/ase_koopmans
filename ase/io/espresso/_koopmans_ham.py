"""Reads Koopmans Ham files

Read structures and results from kcw.x (ham mode) output files.
Read structures from kcw.x (ham mode) input files.
"""

import copy
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.dft.kpoints import BandPath
from ase.utils import basestring
from ._utils import construct_kpoints_card, generic_construct_namelist, safe_string_to_list_of_floats, time_to_float, \
    read_fortran_namelist
from ._wann2kc import KEYS as W2KCW_KEYS

from ase.calculators.espresso import KoopmansHam

KEYS = copy.deepcopy(W2KCW_KEYS)
KEYS['HAM'] = ['do_bands', 'use_ws_distance', 'write_hr', 'l_alpha_corr']


def write_koopmans_ham_in(fd, atoms, input_data=None, pseudopotentials=None,
                          kspacing=None, kpts=None, koffset=(0, 0, 0), **kwargs):

    if 'input_data' in atoms.calc.parameters and input_data is None:
        input_data = atoms.calc.parameters['input_data']

    input_parameters = construct_namelist(input_data, **kwargs)

    assert input_parameters['CONTROL']['calculation'] == 'ham'

    lines = []
    for section in input_parameters:
        assert section in KEYS.keys()

        if section == 'WANNIER' and not input_parameters['CONTROL'].get('kcw_at_ks', True):
            # Do not write the WANNIER section if kcw_at_ks is true
            continue

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

    if input_parameters['HAM'].get('do_bands', True):
        assert isinstance(kpts, BandPath)
        lines += construct_kpoints_card(atoms, kpts, kspacing, koffset)

    fd.writelines(lines)


def read_koopmans_ham_in(fileobj):
    data, _ = read_fortran_namelist(fileobj)
    calc = KoopmansHam(**{k: v for block in data.values() for k, v in block.items()})
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
    walltime = None
    kpts = []
    eigenvalues = []
    ks_eigenvalues_on_grid = []
    ki_eigenvalues_on_grid = []

    for i_line, line in enumerate(flines):

        if 'KC interpolated eigenvalues at k=' in line:
            kpts.append([float(x) for x in line.split()[-3:]])
            eigenvalues.append([])
            j_line = i_line + 2
            while True:
                line2 = flines[j_line].strip()
                if len(line2) == 0:
                    break
                eigenvalues[-1] += safe_string_to_list_of_floats(line2)
                j_line += 1

        if 'INFO: KI[2nd] HAMILTONIAN CALCULATION ik=' in line:
            ks_eigenvalues_on_grid.append([])
            ki_eigenvalues_on_grid.append([])

        if line.startswith('          KS '):
            ks_eigenvalues_on_grid[-1] += safe_string_to_list_of_floats(line.replace('KS', ''))

        if line.startswith('          KI '):
            ki_eigenvalues_on_grid[-1] += safe_string_to_list_of_floats(line.replace('KI', ''))

        if 'JOB DONE' in line:
            job_done = True

        if 'KC_WANN      :' in line:
            time_str = line.split()[-2]
            walltime = time_to_float(time_str)

    # If doing a gamma-point-only calculation, populate the eigenvalues (which won't be printed)
    if len(eigenvalues) == 0 and len(ki_eigenvalues_on_grid) == 1:
        eigenvalues = ki_eigenvalues_on_grid

    # Put everything together
    calc = SinglePointDFTCalculator(structure)
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    calc.results['eigenvalues'] = eigenvalues
    calc.results['ks_eigenvalues_on_grid'] = ks_eigenvalues_on_grid
    calc.results['ki_eigenvalues_on_grid'] = ki_eigenvalues_on_grid
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=False, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
