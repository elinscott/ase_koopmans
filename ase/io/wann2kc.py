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

from ase.io.espresso import Namelist, read_espresso_in, write_espresso_in
from ase.io.espresso import construct_namelist as espresso_construct_namelist

from ase.calculators.wann2kc import Wann2KC

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

KEYS = Namelist((
    ('control', ['prefix', 'outdir', 'kc_iverbosity', 'kc_at_ks', 'homo_only', 'read_unitary_matrix', 'l_vcut']),
    ('wannier', ['seedname', 'check_ks', 'num_wann_occ', 'num_wann_emp', 'have_empty', 'has_disentangle'])))


def write_wann2kc_in(fd, atoms, input_data=None, **kwargs):

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
    '''
    Using espresso.py's construct_namelist function with differently defined "KEYS"
    '''
    return espresso_construct_namelist(parameters, warn, KEYS, **kwargs)
