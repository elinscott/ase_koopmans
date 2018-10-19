"""
Stream input commands to lammps to perform desired simulations
"""
from ase.parallel import paropen
from ase.utils import basestring as asestring
from ase.calculators.lammps import Prism

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = '__end_of_ase_invoked_calculation__'


def write_lammps_in(lammps_in, parameters, atoms, lammps_trj=None,
                    lammps_data=None, verbose=True, prismobj=None):
    """Write a LAMMPS in_ file with run parameters and settings."""

    if isinstance(lammps_in, asestring):
        f = paropen(lammps_in, 'wb')
        close_in_file = True
    else:
        # Expect lammps_in to be a file-like object
        f = lammps_in
        close_in_file = False


    if verbose:
        f.write('# (written by ASE)\n'.encode('utf-8'))

    # Write variables
    f.write(('clear\n'
             'variable dump_file string "{0}"\n'
             'variable data_file string "{1}"\n'
            ).format(lammps_trj, lammps_data).encode('utf-8'))

    if 'package' in parameters:
        f.write(('\n'.join(['package {0}'.format(p)
                            for p in parameters['package']]) +
                 '\n').encode('utf-8'))

    pbc = atoms.get_pbc()
    f.write('units metal \n'.encode('utf-8'))
    if 'boundary' in parameters:
        f.write('boundary {0} \n'.format(
            parameters['boundary']).encode('utf-8'))
    else:
        f.write('boundary {0} {1} {2} \n'.format(
            *tuple('sp'[int(x)] for x in pbc)).encode('utf-8'))
    f.write('atom_modify sort 0 0.0 \n'.encode('utf-8'))
    for key in ('neighbor', 'newton'):
        if key in parameters:
            f.write('{0} {1} \n'.format(
                key, parameters[key]).encode('utf-8'))
    f.write('\n'.encode('utf-8'))

    # write the simulation box and the atoms
    if not lammps_data:
        if verbose:
            f.write('## Original ase cell\n'.encode('utf-8'))
            f.write(''.join(['# {0:.16} {1:.16} {2:.16}\n'.format(*x)
                             for x in atoms.get_cell()]
                           ).encode('utf-8'))

        f.write('lattice sc 1.0\n'.encode('utf-8'))
        if not prismobj:
            prismobj = Prism(atoms.get_cell())
        xhi, yhi, zhi, xy, xz, yz = prismobj.get_lammps_prism_str()
        if parameters['always_triclinic'] or prismobj.is_skewed():
            f.write('region asecell prism 0.0 {0} 0.0 {1} 0.0 {2} '
                    ''.format(xhi, yhi, zhi).encode('utf-8'))
            f.write('{0} {1} {2} side in units box\n'
                    ''.format(xy, xz, yz).encode('utf-8'))
        else:
            f.write(('region asecell block 0.0 {0} 0.0 {1} 0.0 {2} '
                     'side in units box\n').format(
                         xhi, yhi, zhi).encode('utf-8'))

        symbols = atoms.get_chemical_symbols()
        try:
            # By request, specific atom type ordering
            species = parameters['specorder']
        except AttributeError:
            # By default, atom types in alphabetic order
            species = sorted(set(symbols))

        n_atom_types = len(species)
        species_i = dict([(s, i + 1) for i, s in enumerate(species)])

        f.write('create_box {0} asecell\n'.format(
            n_atom_types).encode('utf-8'))
        for s, pos in zip(symbols, atoms.get_positions()):
            if verbose:
                f.write('# atom pos in ase cell: {0:.16} {1:.16} {2:.16}'
                        '\n'.format(*tuple(pos)).encode('utf-8'))
            f.write('create_atoms {0} single {1} {2} {3} units box\n'
                    ''.format(
                        *((species_i[s],) +
                          prismobj.pos_to_lammps_fold_str(pos))
                    ).encode('utf-8'))

    # or simply refer to the data-file
    else:
        f.write('read_data {0}\n'.format(lammps_data).encode('utf-8'))

    # Write interaction stuff
    f.write('\n### interactions \n'.encode('utf-8'))
    if ('pair_style' in parameters) and ('pair_coeff' in parameters):
        pair_style = parameters['pair_style']
        f.write('pair_style {0} \n'.format(pair_style).encode('utf-8'))
        for pair_coeff in parameters['pair_coeff']:
            f.write('pair_coeff {0} \n'
                    ''.format(pair_coeff).encode('utf-8'))
        if 'mass' in parameters:
            for mass in parameters['mass']:
                f.write('mass {0} \n'.format(mass).encode('utf-8'))
    else:
        # simple default parameters
        # that should always make the LAMMPS calculation run
        f.write('pair_style lj/cut 2.5 \n'
                'pair_coeff * * 1 1 \n'
                'mass * 1.0 \n'.encode('utf-8'))

    if 'group' in parameters:
        f.write(('\n'.join(['group {0}'.format(p)
                            for p in parameters['group']]) +
                 '\n').encode('utf-8'))

    f.write(
        '\n### run\n'
        'fix fix_nve all nve\n'.encode('utf-8'))

    if 'fix' in parameters:
        f.write(('\n'.join(['fix {0}'.format(p)
                            for p in parameters['fix']]) +
                 '\n').encode('utf-8'))

    f.write(
        'dump dump_all all custom {1} {0} id type x y z vx vy vz '
        'fx fy fz\n'
        ''.format(lammps_trj, parameters['dump_period']).encode('utf-8'))
    f.write('thermo_style custom {0}\n'
            'thermo_modify flush yes\n'
            'thermo 1\n'.format(
                ' '.join(parameters['thermo_args'])).encode('utf-8'))

    if 'timestep' in parameters:
        f.write('timestep {0}\n'.format(
            parameters['timestep']).encode('utf-8'))

    if 'minimize' in parameters:
        f.write('minimize {0}\n'.format(
            parameters['minimize']).encode('utf-8'))
    if 'run' in parameters:
        f.write('run {0}\n'.format(parameters['run']).encode('utf-8'))
    if not (('minimize' in parameters) or ('run' in parameters)):
        f.write('run 0\n'.encode('utf-8'))

    f.write('print "{0}" \n'.format(CALCULATION_END_MARK).encode('utf-8'))
    # Force LAMMPS to flush log
    f.write('log /dev/stdout\n'.encode('utf-8'))

    f.flush()
    if close_in_file:
        f.close()
