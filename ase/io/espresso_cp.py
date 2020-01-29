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

# For reading in JSON files containing CP directives
import json
from ase.io.jsonio import numpyfy

import numpy as np

from ase.atoms import Atoms
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)
# from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
# from ase.dft.kpoints import kpoint_convert
# from ase.constraints import FixAtoms, FixCartesian
# from ase.data import chemical_symbols, atomic_numbers
from ase.units import create_units, Hartree, Bohr
from ase.utils import basestring

from ase.io.espresso import Namelist, KEYS, SSSP_VALENCE, \
   read_espresso_in, ibrav_to_cell, get_atomic_positions, \
   get_cell_parameters, str_to_value, read_fortran_namelist, ffloat, \
   label_to_symbol, infix_float, construct_namelist, grep_valence, \
   cell_to_ibrav, kspacing_to_grid, write_espresso_in, get_constraint

from ase.calculators.espresso_cp import Espresso_cp

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

KEYS['CONTROL']   += ['ndr', 'ndw']
KEYS['SYSTEM']    += ['fixed_band', 'f_cutoff', 'restart_from_wannier_pwscf', 'do_orbdep', 
                      'fixed_state', 'do_ee', 'nelec', 'nelup', 'neldw']
KEYS['ELECTRONS'] += ['empty_states_nbnd', 'maxiter', 'empty_states_maxstep']
KEYS['NKSIC']      = ['do_innerloop', 'one_innerloop_only', 'nkscalfact', 'odd_nkscalfact', 
                      'odd_nkscalfact_empty', 'which_orbdep', 'print_wfc_anion', 
                      'index_empty_to_save', 'esic_conv_thr']

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

def read_espresso_cp_json(fileobj):
    """Reads JSON files containing CP directives

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename

    Returns
    -------
    A dictionary of atoms objects with Espresso_cp calculators

    """
    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')

    # Decode JSON file
    cp_decoder = json.JSONDecoder().decode
    bigdct = numpyfy(cp_decoder(fileobj.read()))

    all_atoms = {}
    for bigkey, dct in bigdct.items():
        dct = {k.lower() : v for k, v in dct.items()}

        # Create a blank atoms object
        calc = Espresso_cp()

        # Storing CP keywords in atoms.calc.parameters
        cp_keys = [key for block in KEYS.values() for key in block]
        cp_dct = {k : v for k, v in dct.items() if k in cp_keys}
        dct = {k : v for k, v in dct.items() if k not in cp_keys}
        calc.parameters['input_data'] = construct_namelist(cp_dct)

        # Going over keywords that correspond to CP blocks
        if 'atomic_positions' in dct:
            dct['positions'] = np.array([l.split()[1:] for l in dct['atomic_positions'].strip().split('\n')], dtype=float)
            del dct['atomic_positions']
        if 'cell_size' in dct:
            dct['cell'] = np.array([l.split() for l in dct['cell_size'].strip().split('\n')], dtype=float)
            del dct['cell_size']
        if 'atomic_species' in dct:
            calc.parameters['pseudopotentials'] = {l.split()[0] : l.split()[2] for l in dct['atomic_species'].strip().split('\n')}
            del dct['atomic_species']
        # Dealing with units
        units = {'bohr' : Bohr, 'ang' : 1.0}
        if 'unit_cell' in dct:
            dct['cell'] *= units[dct['unit_cell'].lower()]
            del dct['unit_cell']

        # Assuming all remaining keys will be parseable by ASE
        atoms = Atoms(**dct)
        calc.atoms = atoms
        atoms.calc = calc

        all_atoms[bigkey] = atoms

    return all_atoms

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
            odd_energy = float(line.split()[3])*Hartree

        if 'HOMO Eigenvalue (eV)' in line:
            homo_energy = float(cpo_lines[i_line + 2])
    
        if 'LUMO Eigenvalue (eV)' in line:
            lumo_energy = float(cpo_lines[i_line + 2])

        if 'JOB DONE' in line:
            job_done = True

        if 'Orb -- Charge  ---' in line:
            for key in orbital_data:
                orbital_data[key].append([])
            j_line = i_line + 2
            line = cpo_lines[j_line]
            while 'OCC' in line:
                splitline = line.split()
                orbital_data['charge'][-1].append(float(splitline[3]))
                orbital_data['centres'][-1].append([float(x)*Bohr for x in splitline[5:8]])
                orbital_data['spreads'][-1].append(float(splitline[9])*Bohr**2)
                orbital_data['self-Hartree'][-1].append(float(splitline[10]))
                j_line += 1
                line = cpo_lines[j_line]

    # Forces
    # forces = None
    # for force_index in indexes[_CP_FORCE]:
    #     if image_index < force_index < next_index:
    #         # Before QE 5.3 'negative rho' added 2 lines before forces
    #         # Use exact lines to stop before 'non-local' forces
    #         # in high verbosity
    #         if not cpo_lines[force_index + 2].strip():
    #             force_index += 4
    #         else:
    #             force_index += 2
    #         # assume contiguous
    #         forces = [
    #             [float(x) for x in force_line.split()[-3:]] for force_line
    #             in cpo_lines[force_index:force_index + len(structure)]]
    #         forces = np.array(forces) * units['Ry'] / units['Bohr']

    # Stress
    # stress = None
    # for stress_index in indexes[_CP_STRESS]:
    #     if image_index < stress_index < next_index:
    #         sxx, sxy, sxz = cpo_lines[stress_index + 1].split()[:3]
    #         _, syy, syz = cpo_lines[stress_index + 2].split()[:3]
    #         _, _, szz = cpo_lines[stress_index + 3].split()[:3]
    #         stress = np.array([sxx, syy, szz, syz, sxz, sxy], dtype=float)
    #         # sign convention is opposite of ase
    #         stress *= -1 * units['Ry'] / (units['Bohr'] ** 3)

    # Magmoms
    # magmoms = None
    # for magmoms_index in indexes[_CP_MAGMOM]:
    #     if image_index < magmoms_index < next_index:
    #         magmoms = [
    #             float(mag_line.split()[5]) for mag_line
    #             in cpo_lines[magmoms_index + 1:
    #                          magmoms_index + 1 + len(structure)]]

    # Fermi level
    # efermi = None
    # for fermi_index in indexes[_CP_FERMI]:
    #     if image_index < fermi_index < next_index:
    #         efermi = float(cpo_lines[fermi_index].split()[-2])

    # K-points
    # ibzkpts = None
    # weights = None
    # kpoints_warning = "Number of k-points >= 100: " + \
    #                   "set verbosity='high' to print them."

    # for kpts_index in indexes[_CP_KPTS]:
    #     nkpts = int(cpo_lines[kpts_index].split()[4])
    #     kpts_index += 2

    #     if cpo_lines[kpts_index].strip() == kpoints_warning:
    #         continue

    #     # QE prints the k-points in units of 2*pi/alat
    #     # with alat defined as the length of the first
    #     # cell vector
    #     cell = structure.get_cell()
    #     alat = np.linalg.norm(cell[0])
    #     ibzkpts = []
    #     weights = []
    #     for i in range(nkpts):
    #         l = cpo_lines[kpts_index + i].split()
    #         weights.append(float(l[-1]))
    #         coord = np.array([l[-6], l[-5], l[-4].strip('),')],
    #                          dtype=float)
    #         coord *= 2 * np.pi / alat
    #         coord = kpoint_convert(cell, ckpts_kv=coord)
    #         ibzkpts.append(coord)
    #     ibzkpts = np.array(ibzkpts)
    #     weights = np.array(weights)

    # kpts = None
    # kpoints_warning = "Number of k-points >= 100: " + \
    #                   "set verbosity='high' to print the bands."

    # for bands_index in indexes[_CP_BANDS] + indexes[_CP_BANDSTRUCTURE]:
    #     if image_index < bands_index < next_index:
    #         bands_index += 2

    #         if cpo_lines[bands_index].strip() == kpoints_warning:
    #             continue

    #         assert ibzkpts is not None
    #         spin, bands, eigenvalues = 0, [], [[], []]

    #         while True:
    #             l = cpo_lines[bands_index].replace('-', ' -').split()
    #             if len(l) == 0:
    #                 if len(bands) > 0:
    #                     eigenvalues[spin].append(bands)
    #                     bands = []
    #             elif l == ['occupation', 'numbers']:
    #                 bands_index += 3
    #             elif l[0] == 'k' and l[1].startswith('='):
    #                 pass
    #             elif 'SPIN' in l:
    #                 if 'DOWN' in l:
    #                     spin += 1
    #             else:
    #                 try:
    #                     bands.extend(map(float, l))
    #                 except ValueError:
    #                     break
    #             bands_index += 1

    #         if spin == 1:
    #             assert len(eigenvalues[0]) == len(eigenvalues[1])
    #         assert len(eigenvalues[0] + eigenvalues[1]) == len(ibzkpts)

    #         kpts = []
    #         for s in range(spin + 1):
    #             for w, k, e in zip(weights, ibzkpts, eigenvalues[s]):
    #                 kpt = SinglePointKPoint(w, s, k, eps_n=e)
    #                 kpts.append(kpt)

    # Put everything together
    calc = SinglePointDFTCalculator(structure, energy=energy) #,
    #                                 forces=forces, stress=stress,
    #                                 magmoms=magmoms, efermi=efermi,
    #                                 ibzkpts=ibzkpts)
    # calc.kpts = kpts
    calc.results['energy'] = energy
    calc.results['odd_energy'] = odd_energy
    calc.results['homo_energy'] = homo_energy
    calc.results['lumo_energy'] = lumo_energy
    calc.results['eigenvalues'] = eigenvalues
    calc.results['lambda_ii'] = lambda_ii
    calc.results['job_done'] = job_done
    calc.results['orbital_data'] = orbital_data
    structure.set_calculator(calc)

    yield structure
