"""Reads Quantum ESPRESSO files.

Read structures and results from kcp.x output files. Read
structures from kcp.x input files.
"""

import copy
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.utils import basestring
from ._utils import generic_construct_namelist, units, time_to_float, safe_string_to_list_of_floats, write_espresso_in
from ._pw import read_pw_in
from ._pw import KEYS as PW_KEYS
from ase.calculators.espresso import Espresso_kcp


KEYS = copy.deepcopy(PW_KEYS)
KEYS['CONTROL'] += ['ndr', 'ndw', 'ekin_conv_thr', 'write_hr', 'print_real_space_density']
KEYS['SYSTEM'] += ['fixed_band', 'f_cutoff', 'restart_from_wannier_pwscf', 'do_orbdep',
                   'fixed_state', 'do_ee', 'nelec', 'nelup', 'neldw', 'do_wf_cmplx',
                   'nr1b', 'nr2b', 'nr3b']
KEYS['ELECTRONS'] += ['maxiter', 'empty_states_maxstep',
                      'electron_dynamics', 'passop', 'do_outerloop', 'do_outerloop_empty']
KEYS['EE'] = ['which_compensation', 'tcc_odd']
KEYS['NKSIC'] = ['do_innerloop', 'nkscalfact', 'odd_nkscalfact',
                 'odd_nkscalfact_empty', 'which_orbdep', 'print_wfc_anion',
                 'index_empty_to_save', 'innerloop_cg_nreset', 'innerloop_cg_nsd',
                 'innerloop_init_n', 'hartree_only_sic', 'esic_conv_thr',
                 'do_innerloop_cg', 'innerloop_nmax', 'do_innerloop_empty',
                 'innerloop_cg_ratio', 'fref', 'kfact', 'wo_odd_in_empty_run',
                 'aux_empty_nbnd', 'print_evc0_occ_empty', 'do_bare_eigs']
KEYS['IONS'] += ['ion_nstepe', 'ion_radius']

# Section identifiers
_CP_START = 'CP: variable-cell Car-Parrinello molecular dynamics'
_CP_END = 'This run was terminated on:'
_CP_CELL = 'CELL_PARAMETERS'
_CP_POS = 'ATOMIC_POSITIONS'
# _CP_MAGMOM =
# _CP_FORCE =
_CP_TOTEN = '                total energy'
_CP_BANDS = 'Eigenvalues (eV), kp'
_CP_LAMBDA = 'fixed_lambda'
# _CP_STRESS =
# _CP_FERMI =
# _CP_KPTS =
# _CP_BANDSTRUCTURE =


def write_koopmans_cp_in(*args, **kwargs):
    kwargs['local_construct_namelist'] = construct_namelist
    kwargs['include_kpoints'] = False
    write_espresso_in(*args, **kwargs)


def read_koopmans_cp_in(fileobj):
    atoms = read_pw_in(fileobj)

    # Generating Espresso_kcp calculator from Espresso calculator
    calc = Espresso_kcp(**atoms.calc.parameters)

    # Overwriting the Espresso calculator with the new Espresso_kcp calculator
    atoms.calc = calc
    atoms.calc.atoms = atoms

    return atoms


def read_koopmans_cp_out(fileobj, index=-1, results_required=True):
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
    mp1_energy = None
    mp2_energy = None
    lambda_ii = None
    eigenvalues = []
    job_done = False
    orbital_data = {'charge': [], 'centres': [], 'spreads': [], 'self-Hartree': []}
    walltime = None
    convergence = {'filled': [], 'empty': []}

    convergence_key = 'filled'
    i_spin_orbital_data = None

    for i_line, line in enumerate(cpo_lines):

        # Energy
        if _CP_TOTEN in line:
            energy = float(line.split()[-3]) * units.Hartree

        if _CP_LAMBDA in line and lambda_ii is None:
            lambda_ii = float(line.split()[-1]) * units.Hartree

        # Bands
        if _CP_BANDS in line:
            if 'Empty' not in line:
                eigenvalues.append([])

            j_line = i_line + 2
            while len(cpo_lines[j_line].strip()) > 0:
                eigenvalues[-1] += safe_string_to_list_of_floats(cpo_lines[j_line])
                j_line += 1

        if 'odd energy' in line:
            odd_energy = float(line.split()[3]) * units.Hartree

        if 'HOMO Eigenvalue (eV)' in line and '*' not in cpo_lines[i_line + 2]:
            homo_energy = float(cpo_lines[i_line + 2])

        if 'LUMO Eigenvalue (eV)' in line and '*' not in cpo_lines[i_line + 2]:
            lumo_energy = float(cpo_lines[i_line + 2])

        if 'Makov-Payne 1-order energy' in line:
            mp1_energy = float(line.split()[4]) * units.Hartree

        if 'Makov-Payne 2-order energy' in line:
            mp2_energy = float(line.split()[4]) * units.Hartree

        if 'JOB DONE' in line:
            job_done = True

        # Start of block of orbitals for a given spin channel
        if 'Orb -- Charge  ---' in line or 'Orb -- Empty Charge' in line:
            i_spin_orbital_data = int(line.split()[-1]) - 1
            for key in orbital_data:
                if len(orbital_data[key]) < i_spin_orbital_data + 1:
                    orbital_data[key].append([])

        # Orbital information
        if line.startswith(('OCC', 'EMP')):
            if 'NaN' in line:
                continue
            line = line.replace('********', '   0.000')
            values = [float(line[i - 4:i + 4]) for i, c in enumerate(line) if c == '.']
            orbital_data['charge'][i_spin_orbital_data].append(values[0])
            orbital_data['centres'][i_spin_orbital_data].append([x * units.Bohr for x in values[1:4]])
            orbital_data['spreads'][i_spin_orbital_data].append(values[4] * units.Bohr**2)
            orbital_data['self-Hartree'][i_spin_orbital_data].append(values[5])

        # Tracking convergence
        if 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EMPTY STATES' in line:
            convergence_key = 'empty'

        if 'iteration = ' in line and 'eff iteration = ' in line:
            values = [l.split()[0] for l in line.split('=')[1:]]
            [it, eff_it, etot] = values[:3]
            entry = {'iteration': int(it), 'eff iteration': int(eff_it),
                     'Etot': float(etot) * units.Hartree}
            if len(values) == 4:
                entry['delta_E'] = float(values[3]) * units.Hartree
            convergence[convergence_key].append(entry)

        if 'wall time' in line:
            time_str = line.split(',')[1].strip().rstrip('wall time')
            walltime = time_to_float(time_str)

    # Put everything together
    calc = SinglePointDFTCalculator(structure, energy=energy)  # ,
    #                                 forces=forces, stress=stress,
    #                                 magmoms=magmoms, efermi=efermi,
    #                                 ibzkpts=ibzkpts)
    calc.results['energy'] = energy
    calc.results['odd_energy'] = odd_energy
    calc.results['homo_energy'] = homo_energy
    calc.results['lumo_energy'] = lumo_energy
    calc.results['mp1_energy'] = mp1_energy
    calc.results['mp2_energy'] = mp2_energy
    calc.results['eigenvalues'] = eigenvalues
    calc.results['lambda_ii'] = lambda_ii
    calc.results['orbital_data'] = orbital_data
    calc.results['convergence'] = convergence
    calc.results['job_done'] = job_done
    calc.results['walltime'] = walltime
    structure.calc = calc

    yield structure


def construct_namelist(parameters=None, warn=False, **kwargs):
    return generic_construct_namelist(parameters, warn, KEYS, **kwargs)
