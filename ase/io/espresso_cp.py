"""Reads Quantum ESPRESSO files.

Read multiple structures and results from pw.x output files. Read
structures from cp.x input files.

Built for CP

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
   get_cell_parameters, str_to_value, read_fortran_namelist, ffloat, \
   label_to_symbol, infix_float, grep_valence, \
   cell_to_ibrav, kspacing_to_grid, write_espresso_in, get_constraint
from ase.io.espresso import KEYS as PW_KEYS

from ase.calculators.espresso_cp import Espresso_cp

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

KEYS = copy.deepcopy(PW_KEYS)
KEYS['CONTROL']   += ['ndr', 'ndw', 'ekin_conv_thr', 'write_hr']
KEYS['SYSTEM']    += ['fixed_band', 'f_cutoff', 'restart_from_wannier_pwscf', 'do_orbdep', 
                      'fixed_state', 'do_ee', 'nelec', 'nelup', 'neldw', 'do_wf_cmplx', 
                      'nr1b', 'nr2b', 'nr3b']
KEYS['ELECTRONS'] += ['empty_states_nbnd', 'maxiter', 'empty_states_maxstep', 
                      'electron_dynamics', 'passop', 'do_outerloop', 'do_outerloop_empty']
KEYS['EE']         = ['which_compensation', 'tcc_odd']
KEYS['NKSIC']      = ['do_innerloop', 'nkscalfact', 'odd_nkscalfact', 
                      'odd_nkscalfact_empty', 'which_orbdep', 'print_wfc_anion', 
                      'index_empty_to_save', 'innerloop_cg_nreset', 'innerloop_cg_nsd', 
                      'innerloop_init_n', 'hartree_only_sic', 'esic_conv_thr', 
                      'do_innerloop_cg', 'innerloop_nmax', 'do_innerloop_empty', 
                      'innerloop_cg_ratio', 'fref', 'kfact', 'wo_odd_in_empty_run', 
                      'aux_empty_nbnd', 'print_evc0_occ_empty']
KEYS['IONS']      += ['ion_nstepe', 'ion_radius(1)', 'ion_radius(2)', 'ion_radius(3)',
                      'ion_radius(4)'] 

# Section identifiers
_CP_START = 'CP: variable-cell Car-Parrinello molecular dynamics'
_CP_END = 'This run was terminated on:'
_CP_CELL = 'CELL_PARAMETERS'
_CP_POS = 'ATOMIC_POSITIONS'
# _CP_MAGMOM =
# _CP_FORCE =
_CP_TOTEN = '                total energy'
_CP_BANDS = '   Eigenvalues (eV), kp'
_CP_LAMBDA = 'fixed_lambda'
# _CP_STRESS =
# _CP_FERMI =
# _CP_KPTS =
# _CP_BANDSTRUCTURE = 

def write_espresso_cp_in(fd, atoms, input_data=None, pseudopotentials=None,
                      kspacing=None, kpts=None, koffset=(0, 0, 0),
                      **kwargs):

    write_espresso_in(fd, atoms, input_data, pseudopotentials,
                      kspacing, kpts, koffset, **kwargs)

    if not fd.closed:
        fd.close()
    
    # Extra blocks
    extra_lines = []
    input_parameters = construct_namelist(input_data, **kwargs)
    for section in input_parameters:
        if section.lower() not in ['ee', 'nksic']:
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
    extra_lines.append('\n')

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


def read_espresso_cp_in(fileobj):
    atoms = read_espresso_in(fileobj)

    # Generating Espresso_cp calculator from Espresso calculator
    data = atoms.calc.parameters['input_data']
    pseudos = atoms.calc.parameters['pseudopotentials']
    calc = Espresso_cp(input_data=data, pseudopotentials=pseudos)

    # Overwriting the Espresso calculator with the new Espresso_cp calculator
    atoms.set_calculator(calc)
    atoms.calc.atoms = atoms

    return atoms


def read_espresso_cp_out(fileobj, index=-1, results_required=True):
    """Reads Quantum ESPRESSO output files.

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
    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')

    # work with a copy in memory for faster random access
    cpo_lines = fileobj.readlines()

    # For the moment, provide an empty atoms object
    structure = Atoms()

    # Extract calculation results
    energy = None
    odd_energy = None
    lumo_energy = None
    homo_energy = None
    lambda_ii = None
    eigenvalues = []
    job_done = False
    orbital_data = {'charge' : [], 'centres' : [], 'spreads' : [], 'self-Hartree' : []}
    walltime = None
    convergence = {'filled' : [], 'empty': []}

    convergence_key = 'filled'

    for i_line, line in enumerate(cpo_lines):

        # Energy
        if _CP_TOTEN in line:
            energy = float(line.split()[-3])*units.Hartree

        if _CP_LAMBDA in line and lambda_ii is None:
            lambda_ii = float(line.split()[-1])*units.Hartree

        # Bands
        if _CP_BANDS in line:
            try:
                eigenvalues.append([float(e) for e in cpo_lines[i_line + 2].split()])
            except:
                pass
            if 'Empty States Eigenvalues' in cpo_lines[i_line + 4]:
                try:
                    eigenvalues[-1] += [float(e) for e in cpo_lines[i_line + 6].split()]
                except:
                    pass

        if 'odd energy' in line:
            odd_energy = float(line.split()[3])*units.Hartree

        if 'HOMO Eigenvalue (eV)' in line and '*' not in cpo_lines[i_line + 2]:
            homo_energy = float(cpo_lines[i_line + 2])
    
        if 'LUMO Eigenvalue (eV)' in line and '*' not in cpo_lines[i_line + 2]:
            lumo_energy = float(cpo_lines[i_line + 2])

        if 'JOB DONE' in line:
            job_done = True

        # Start of block of orbitals for a given spin channel
        if 'Orb -- Charge  ---' in line:
            for key in orbital_data:
                orbital_data[key].append([])

        # Orbital information
        if line.startswith(('OCC', 'EMP')):
            if 'NaN' in line:
                continue
            line = line.replace('********', '   0.000')
            values = [float(line[i-4:i+4]) for i, c in enumerate(line) if c == '.']
            orbital_data['charge'][-1].append(values[0])
            orbital_data['centres'][-1].append([x*units.Bohr for x in values[1:4]])
            orbital_data['spreads'][-1].append(values[4]*units.Bohr**2)
            orbital_data['self-Hartree'][-1].append(values[5])

        # Tracking convergence
        if 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EMPTY STATES' in line:
            convergence_key = 'empty'

        if 'iteration = ' in line and 'eff iteration = ' in line:
            values = [l.split()[0] for l in line.split('=')[1:]]
            [it, eff_it, etot] = values[:3]
            entry = {'iteration': int(it), 'eff iteration': int(eff_it), 
                     'Etot': float(etot)*units.Hartree}
            if len(values) == 4:
                entry['delta_E'] = float(values[3])*units.Hartree
            convergence[convergence_key].append(entry)

        if 'wall time' in line:
            time_str = line.split(',')[1].strip().rstrip('wall time')
            if 'h' in time_str:
                hours, rem = time_str.split('h')
            else:
                hours, rem = 0, time_str
            if 'm' in rem:
                minutes, rem = rem.split('m')
            else:
                minutes = 0
            if 's' in rem:
                seconds = rem.rstrip('s')
            else:
                seconds = 0
            walltime = (float(hours)*60 + float(minutes))*60 + float(seconds)

    # Put everything together
    calc = SinglePointDFTCalculator(structure, energy=energy) #,
    #                                 forces=forces, stress=stress,
    #                                 magmoms=magmoms, efermi=efermi,
    #                                 ibzkpts=ibzkpts)
    calc.results['energy'] = energy
    calc.results['odd_energy'] = odd_energy
    calc.results['homo_energy'] = homo_energy
    calc.results['lumo_energy'] = lumo_energy
    calc.results['eigenvalues'] = eigenvalues
    calc.results['lambda_ii'] = lambda_ii
    calc.results['orbital_data'] = orbital_data
    calc.results['convergence'] = convergence
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    structure.set_calculator(calc)

    yield structure


def construct_namelist(parameters=None, warn=False, **kwargs):
    """
    Construct an ordered Namelist containing all the parameters given (as
    a dictionary or kwargs). Keys will be inserted into their appropriate
    section in the namelist and the dictionary may contain flat and nested
    structures. Any kwargs that match input keys will be incorporated into
    their correct section. All matches are case-insensitive, and returned
    Namelist object is a case-insensitive dict.

    If a key is not known to ase, but in a section within `parameters`,
    it will be assumed that it was put there on purpose and included
    in the output namelist. Anything not in a section will be ignored (set
    `warn` to True to see ignored keys).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    The priority of the keys is:
        kwargs[key] > parameters[key] > parameters[section][key]
    Only the highest priority item will be included.

    Copied from ase/io/espresso.cp

    Parameters
    ----------
    parameters: dict
        Flat or nested set of input parameters.
    warn: bool
        Enable warnings for unused keys.

    Returns
    -------
    input_namelist: Namelist
        cp.x compatible namelist of input parameters.

    """
    # Convert everything to Namelist early to make case-insensitive
    if parameters is None:
        parameters = Namelist()
    else:
        # Maximum one level of nested dict
        # Don't modify in place
        parameters_namelist = Namelist()
        for key, value in parameters.items():
            if isinstance(value, dict):
                parameters_namelist[key] = Namelist(value)
            else:
                parameters_namelist[key] = value
        parameters = parameters_namelist

    # Just a dict
    kwargs = Namelist(kwargs)

    # Final parameter set
    input_namelist = Namelist()

    # Collect
    for section in KEYS:
        sec_list = Namelist()
        for key in KEYS[section]:
            # Check all three separately and pop them all so that
            # we can check for missing values later
            if key in parameters.get(section, {}):
                sec_list[key] = parameters[section].pop(key)
            if key in parameters:
                sec_list[key] = parameters.pop(key)
            if key in kwargs:
                sec_list[key] = kwargs.pop(key)

            # Check if there is a key(i) version (no extra parsing)
            cp_parameters = parameters.copy()
            for arg_key in cp_parameters:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = parameters.pop(arg_key)
            cp_kwargs = kwargs.copy()
            for arg_key in cp_kwargs:
                if arg_key.split('(')[0].strip().lower() == key.lower():
                    sec_list[arg_key] = kwargs.pop(arg_key)

        # Add to output
        input_namelist[section] = sec_list

    unused_keys = list(kwargs)
    # pass anything else already in a section
    for key, value in parameters.items():
        if key in KEYS and isinstance(value, dict):
            input_namelist[key].update(value)
        elif isinstance(value, dict):
            unused_keys.extend(list(value))
        else:
            unused_keys.append(key)

    if warn and unused_keys:
        warnings.warn('Unused keys: {}'.format(', '.join(unused_keys)))

    return input_namelist
