"""Reads Koopmans Ham files

Read structures and results from koopmans_ham.x output files. Read
structures from koopmans_ham.x input files.

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

from ase.io.espresso import read_espresso_in, write_espresso_in, Namelist
from ase.io.espresso import construct_namelist as espresso_construct_namelist
from ase.io.wann2kc import KEYS as W2KKEYS

from ase.calculators.koopmans_ham import KoopmansHam

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

KEYS = copy.deepcopy(W2KKEYS)
KEYS['HAM'] = ['do_bands', 'use_ws_distance', 'write_hr', 'l_alpha_corr', 'lrpa', 'mp1', 'mp2', 'mp3']


def write_koopmans_ham_in(fd, atoms, input_data=None, pseudopotentials=None,
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


def read_koopmans_ham_in(fileobj):
    atoms = read_espresso_in(fileobj)

    # Generating KoopmansHam calculator from Espresso calculator
    data = atoms.calc.parameters['input_data']
    pseudos = atoms.calc.parameters['pseudopotentials']
    calc = KoopmansHam(input_data=data, pseudopotentials=pseudos)

    # Overwriting the Espresso calculator with the new KoopmansHam calculator
    atoms.calc = calc
    atoms.calc.atoms = atoms

    raise NotImplementedError('Yet to test this function')

    return atoms


def read_koopmans_ham_out(fileobj, index=-1, results_required=True):
    """Reads Koopmans Ham output files.

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


def construct_namelist(parameters=None, warn=False, **kwargs):
    '''
    Using espresso.py's construct_namelist function with differently defined "KEYS"
    '''
    espresso_construct_namelist(parameters, warn, KEYS, **kwargs)
