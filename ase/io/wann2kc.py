"""Reads Wannier to KC files

Read structures and results from koopmans_ham.x output files. Read
structures from wannier_to_kc.x input files.

Units are converted using CODATA 2006, as used internally by Quantum
ESPRESSO.
"""

import os
import operator as op
import warnings
from collections import OrderedDict
from os import path
import copy

import numpy as np

from ase.atoms import Atoms
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)
# from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
# from ase.dft.kpoints import kpoint_convert
# from ase.constraints import FixAtoms, FixCartesian
# from ase.data import chemical_symbols, atomic_numbers
from ase.units import create_units
from ase.utils import basestring

from ase.io.espresso import Namelist, SSSP_VALENCE, \
    read_espresso_in, ibrav_to_cell, get_atomic_positions, \
    get_cell_parameters, str_to_value, ffloat, \
    label_to_symbol, infix_float, grep_valence, \
    cell_to_ibrav, kspacing_to_grid, write_espresso_in, get_constraint

from ase.calculators.wann2kc import Wann2KC

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

KEYS = Namelist((
    ('control', ['prefix', 'outdir', 'kc_iverbosity ', 'kc_at_ks', 'homo_only ', 'read_unitary_matrix ', 'l_vcut ']),
    ('wannier', ['seedname', 'check_ks', 'num_wann_occ', 'num_wann_emp', 'have_empty', 'has_disentangle'])))


def write_wann2kc_in(fd, atoms, input_data=None, pseudopotentials=None,
                     kspacing=None, kpts=None, koffset=(0, 0, 0), **kwargs):

    raise NotImplementedError('Yet to complete this function')

    write_espresso_in(fd, atoms, input_data, pseudopotentials,
                      kspacing, kpts, koffset, **kwargs)

    if not fd.closed:
        fd.close()

    # Extra blocks
    extra_lines = []
    input_parameters = construct_namelist(input_data, **kwargs)
    for section in input_parameters:
        if section.lower() not in []:
            continue
        extra_lines.append('&{0}\n'.format(section.upper()))
        for key, value in input_parameters[section].items():
            if value is True:
                extra_lines.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                extra_lines.append('   {0:16} = .false.\n'.format(key))
            elif value is not None:
                # repr format to get quotes around strings
                extra_lines.append('   {0:16} = {1!r:}\n'.format(key, value))
        extra_lines.append('/\n')  # terminate section

    # Read in the original file without the NKSIC and EE blocks
    with open(fd.name, 'r') as fd_read:
        lines = fd_read.readlines()

    # Find where to insert the extra blocks
    i_break = lines.index('\n')
    before = lines[:i_break]
    after = lines[i_break:]

    # Rewrite the file with the extra blocks
    with open(fd.name, 'w') as fd_rewrite:
        fd_rewrite.writelines(before + extra_lines + after)


def read_wann2kc_in(fileobj):
    atoms = read_espresso_in(fileobj)

    # Generating Wann2KC calculator from Espresso calculator
    data = atoms.calc.parameters['input_data']
    pseudos = atoms.calc.parameters['pseudopotentials']
    calc = Wann2KC(input_data=data, pseudopotentials=pseudos)

    # Overwriting the Espresso calculator with the new Wann2KC calculator
    atoms.calc = calc
    atoms.calc.atoms = atoms

    raise NotImplementedError('Yet to test this function')

    return atoms


def read_wann2kc_out(fileobj, index=-1, results_required=True):
    """Reads Wannier to KC output files.

    The atomistic configurations as well as results (energy, force, stress,
    magnetic moments) of the calculation are read for all configurations
    within the output file.

    Will probably raise errors for broken or incomplete files.

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename
    index : slice
        The index of configurations to extract.
    results_required : bool
        If True, atomistic configurations that do not have any
        associated results will not be included. This prevents double
        printed configurations and incomplete calculations from being
        returned as the final configuration with no results data.

    Yields
    ------
    structure : Atoms
        The next structure from the index slice. The Atoms has a
        SinglePointCalculator attached with any results parsed from
        the file.

    """

    raise NotImplementedError('Yet to write this function')
    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')

    # work with a copy in memory for faster random access
    flines = fileobj.readlines()

    # For the moment, provide an empty atoms object
    structure = Atoms()

    # Extract calculation results
    energy = None
    job_done = False
    for i_line, line in enumerate(flines):
        continue

    # Put everything together
    calc = SinglePointDFTCalculator(structure, energy=energy)  # ,
    #                                 forces=forces, stress=stress,
    #                                 magmoms=magmoms, efermi=efermi,
    #                                 ibzkpts=ibzkpts)
    calc.results['energy'] = energy
    calc.results['job_done'] = job_done
    structure.calc = calc

    yield structure


# def construct_namelist(parameters=None, warn=False, **kwargs):
#     """
#     Construct an ordered Namelist containing all the parameters given (as
#     a dictionary or kwargs). Keys will be inserted into their appropriate
#     section in the namelist and the dictionary may contain flat and nested
#     structures. Any kwargs that match input keys will be incorporated into
#     their correct section. All matches are case-insensitive, and returned
#     Namelist object is a case-insensitive dict.
#
#     If a key is not known to ase, but in a section within `parameters`,
#     it will be assumed that it was put there on purpose and included
#     in the output namelist. Anything not in a section will be ignored (set
#     `warn` to True to see ignored keys).
#
#     Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
#     so the `i` should be made to match the output.
#
#     The priority of the keys is:
#         kwargs[key] > parameters[key] > parameters[section][key]
#     Only the highest priority item will be included.
#
#     Copied from ase/io/espresso
#
#     Parameters
#     ----------
#     parameters: dict
#         Flat or nested set of input parameters.
#     warn: bool
#         Enable warnings for unused keys.
#
#     Returns
#     -------
#     input_namelist: Namelist
#         wannier_to_kc.x compatible namelist of input parameters.
#
#     """
#     # Convert everything to Namelist early to make case-insensitive
#     if parameters is None:
#         parameters = Namelist()
#     else:
#         # Maximum one level of nested dict
#         # Don't modify in place
#         parameters_namelist = Namelist()
#         for key, value in parameters.items():
#             if isinstance(value, dict):
#                 parameters_namelist[key] = Namelist(value)
#             else:
#                 parameters_namelist[key] = value
#         parameters = parameters_namelist
#
#     # Just a dict
#     kwargs = Namelist(kwargs)
#
#     # Final parameter set
#     input_namelist = Namelist()
#
#     # Collect
#     for section in KEYS:
#         sec_list = Namelist()
#         for key in KEYS[section]:
#             # Check all three separately and pop them all so that
#             # we can check for missing values later
#             if key in parameters.get(section, {}):
#                 sec_list[key] = parameters[section].pop(key)
#             if key in parameters:
#                 sec_list[key] = parameters.pop(key)
#             if key in kwargs:
#                 sec_list[key] = kwargs.pop(key)
#
#             # Check if there is a key(i) version (no extra parsing)
#             kcp_parameters = parameters.copy()
#             for arg_key in kcp_parameters:
#                 if arg_key.split('(')[0].strip().lower() == key.lower():
#                     sec_list[arg_key] = parameters.pop(arg_key)
#             kcp_kwargs = kwargs.copy()
#             for arg_key in kcp_kwargs:
#                 if arg_key.split('(')[0].strip().lower() == key.lower():
#                     sec_list[arg_key] = kwargs.pop(arg_key)
#
#         # Add to output
#         input_namelist[section] = sec_list
#
#     unused_keys = list(kwargs)
#     # pass anything else already in a section
#     for key, value in parameters.items():
#         if key in KEYS and isinstance(value, dict):
#             input_namelist[key].update(value)
#         elif isinstance(value, dict):
#             unused_keys.extend(list(value))
#         else:
#             unused_keys.append(key)
#
#     if warn and unused_keys:
#         warnings.warn('Unused keys: {}'.format(', '.join(unused_keys)))
#
#     return input_namelist
