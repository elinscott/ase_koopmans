"""ASE LAMMPS Calculator Library Version"""


import ctypes
import operator

import numpy as np
from numpy.linalg import norm

from ase_koopmans.calculators.calculator import Calculator
from ase_koopmans.data import (atomic_numbers as ase_koopmans_atomic_numbers,
                      chemical_symbols as ase_koopmans_chemical_symbols,
                      atomic_masses as ase_koopmans_atomic_masses)
from ase_koopmans.calculators.lammps import convert
from ase_koopmans.geometry import wrap_positions

# TODO
# 1. should we make a new lammps object each time ?
# 4. need a routine to get the model back from lammps
# 5. if we send a command to lmps directly then the calculator does
#    not know about it and the energy could be wrong.
# 6. do we need a subroutine generator that converts a lammps string
#   into a python function that can be called
# 8. make matscipy as fallback
# 9. keep_alive not needed with no system changes
#10. it may be a good idea to unify the cell handling with the one found in
#    lammpsrun.py


# this one may be moved to some more generic place
def is_upper_triangular(arr, atol=1e-8):
    """test for upper triangular matrix base_koopmansd on numpy"""
    # must be (n x n) matrix
    assert len(arr.shape) == 2
    assert arr.shape[0] == arr.shape[1]
    return np.allclose(np.tril(arr, k=-1), 0., atol=atol) and \
        np.all(np.diag(arr) >= 0.0)


def convert_cell(ase_koopmans_cell):
    """
    Convert a parallelepiped (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the base_koopmanss are column vectors
    """
    cell = ase_koopmans_cell.T

    if not is_upper_triangular(cell):
        # rotate base_koopmanss into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = cell[:, 0]
        B = cell[:, 1]
        C = cell[:, 2]
        tri_mat[0, 0] = norm(A)
        Ahat = A / norm(A)
        AxBhat = np.cross(A, B) / norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(ase_koopmans_cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans /= volume
        coord_transform = np.dot(tri_mat, trans)

        return tri_mat, coord_transform
    else:
        return cell, None


class LAMMPSlib(Calculator):
    r"""
**Introduction**

LAMMPSlib is an interface and calculator for LAMMPS_. LAMMPSlib uses
the python interface that comes with LAMMPS to solve an atoms model
for energy, atom forces and cell stress. This calculator creates a
'.lmp' object which is a running lammps program, so further commands
can be sent to this object executed until it is explicitly closed. Any
additional variables calculated by lammps can also be extracted. This
is still experimental code.

**Arguments**

====================  ==========================================================
Keyword                                  Description
====================  ==========================================================
``lmpcmds``           list of strings of LAMMPS commands. You need to supply
                      enough to define the potential to be used e.g.

                      ["pair_style eam/alloy",
                       "pair_coeff * * potentials/NiAlH_jea.eam.alloy Ni Al"]

``atom_types``        dictionary of ``atomic_symbol :lammps_atom_type`` pairs,
                      e.g. ``{'Cu':1}`` to bind copper to lammps atom type 1.
                      If <None>, autocreated by assigning lammps atom types in
                      order that they appear in the first used atoms object.

``atom_type_masses``  dictionary of ``atomic_symbol :mass`` pairs,
                      e.g. ``{'Cu':63.546}`` to optionally assign masses that
                      override default ase_koopmans.data.atomic_masses.  Note that since
                      unit conversion is done automatically in this module,
                      these quantities must be given in the standard ase_koopmans mass
                      units (g/mol)

``log_file``          string
                      path to the desired LAMMPS log file

``lammps_header``     string to use for lammps setup. Default is to use
                      metal units and simple atom simulation.

                      lammps_header=['units metal',
                                  'atom_style atomic',
                                  'atom_modify map array sort 0 0'])

``amendments``        extra list of strings of LAMMPS commands to be run post
                      post initialization. (Use: Initialization amendments) e.g.

                      ["mass 1 58.6934"]

``keep_alive``        Boolean
                      whether to keep the lammps routine alive for more commands

====================  ==========================================================


**Requirements**

To run this calculator you must have LAMMPS installed and compiled to
enable the python interface. See the LAMMPS manual.

If the following code runs then lammps is installed correctly.

   >>> from lammps import lammps
   >>> lmp = lammps()

The version of LAMMPS is also important. LAMMPSlib is suitable for
versions after approximately 2011. Prior to this the python interface
is slightly different from that used by LAMMPSlib. It is not difficult
to change to the earlier format.

**LAMMPS and LAMMPSlib**

The LAMMPS calculator is another calculator that uses LAMMPS (the
program) to calculate the energy by generating input files and running
a separate LAMMPS job to perform the analysis. The output data is then
read back into python. LAMMPSlib makes direct use of the LAMMPS (the
program) python interface. As well as directly running any LAMMPS
command line it allows the values of any of LAMMPS variables to be
extracted and returned to python.

**Example**

Provided that the respective potential file is in the working directory, one
can simply run (note that LAMMPS needs to be compiled to work with EAM
potentials)

::

    from ase_koopmans import Atom, Atoms
    from ase_koopmans.build import bulk
    from ase_koopmans.calculators.lammpslib import LAMMPSlib

    cmds = ["pair_style eam/alloy",
            "pair_coeff * * NiAlH_jea.eam.alloy Ni H"]

    Ni = bulk('Ni', cubic=True)
    H = Atom('H', position=Ni.cell.diagonal()/2)
    NiH = Ni + H

    lammps = LAMMPSlib(lmpcmds=cmds, log_file='test.log')

    NiH.calc = lammps
    print("Energy ", NiH.get_potential_energy())


**Implementation**

LAMMPS provides a set of python functions to allow execution of the
underlying C++ LAMMPS code. The functions used by the LAMMPSlib
interface are::

    from lammps import lammps

    lmp = lammps(cmd_args) # initiate LAMMPS object with command line args

    lmp.scatter_atoms('x',1,3,positions) # atom coords to LAMMPS C array
    lmp.command(cmd) # executes a one line cmd string
    lmp.extract_variable(...) # extracts a per atom variable
    lmp.extract_global(...) # extracts a global variable
    lmp.close() # close the lammps object

For a single Ni atom model the following lammps file commands would be run
by invoking the get_potential_energy() method::

    units metal
    atom_style atomic
    atom_modify map array sort 0 0

    region cell prism 0 xhi 0 yhi 0 zhi xy xz yz units box
    create_box 1 cell
    create_atoms 1 single 0 0 0 units box
    mass * 1.0

    ## user lmpcmds get executed here
    pair_style eam/alloy
    pair_coeff * * NiAlH_jea.eam.alloy Ni
    ## end of user lmmpcmds

    run 0

where xhi, yhi and zhi are the lattice vector lengths and xy,
xz and yz are the tilt of the lattice vectors, all to be edited.


**Notes**

.. _LAMMPS: http://lammps.sandia.gov/

* Units: The default lammps_header sets the units to Angstrom and eV
  and for compatibility with ASE Stress is in GPa.

* The global energy is currently extracted from LAMMPS using
  extract_variable since lammps.lammps currently extract_global only
  accepts the following ['dt', 'boxxlo', 'boxxhi', 'boxylo', 'boxyhi',
  'boxzlo', 'boxzhi', 'natoms', 'nlocal'].

* If an error occurs while lammps is in control it will crash
  Python. Check the output of the log file to find the lammps error.

* If the are commands directly sent to the LAMMPS object this may
  change the energy value of the model. However the calculator will not
  know of it and still return the original energy value.

"""

    implemented_properties = ['energy', 'forces', 'stress']

    started = False
    initialized = False

    default_parameters = dict(
        atom_types=None,
        atom_type_masses=None,
        log_file=None,
        lammps_name='',
        keep_alive=False,
        lammps_header=['units metal',
                       'atom_style atomic',
                       'atom_modify map array sort 0 0'],
        amendments=None,
        boundary=True,
        create_box=True,
        create_atoms=True,
        read_molecular_info=False,
        comm=None)

    def __init__(self, *args, **kwargs):
        Calculator.__init__(self, *args, **kwargs)
        self.lmp = None

    def __del__(self):
        if self.started:
            self.lmp.close()
            self.started = False
            self.lmp = None

    def set_cell(self, atoms, change=False):
        lammps_cell, self.coord_transform = convert_cell(atoms.get_cell())

        xhi, xy, xz, _, yhi, yz, _, _, zhi = convert(
                lammps_cell.flatten(order='C'), "distance", "ASE", self.units)
        box_hi = [xhi, yhi, zhi]

        if change:
            cell_cmd = ('change_box all     '
                        'x final 0 {} y final 0 {} z final 0 {}      '
                        'xy final {} xz final {} yz final {} units box'
                        ''.format(xhi, yhi, zhi, xy, xz, yz))
        else:
            # just in case_koopmans we'll want to run with a funny shape box,
            # and here command will only happen once, and before
            # any calculation
            if self.parameters.create_box:
                self.lmp.command('box tilt large')

            # Check if there are any indefinite boundaries. If so, shrink-wrapping will
            # end up being used, but we want to define the LAMMPS region and box fairly
            # tight around the atoms to avoid losing any
            lammps_boundary_conditions = self.lammpsbc(atoms).split()
            if 's' in lammps_boundary_conditions:
                pos = atoms.get_positions()
                if self.coord_transform is not None:
                    pos = np.dot(self.coord_transform, pos.transpose())
                    pos = pos.transpose()
                posmin = np.amin(pos, axis=0)
                posmax = np.amax(pos, axis=0)

                for i in range(0,3):
                    if lammps_boundary_conditions[i] == 's':
                        box_hi[i] = 1.05*abs(posmax[i] - posmin[i])

            cell_cmd = ('region cell prism    '
                        '0 {} 0 {} 0 {}     '
                        '{} {} {}     units box'
                        ''.format(*box_hi, xy, xz, yz))

        self.lmp.command(cell_cmd)

    def set_lammps_pos(self, atoms):
        # Create local copy of positions that are wrapped along any periodic
        # directions
        pos = wrap_positions(atoms.get_positions(), atoms.get_cell(),
                             atoms.get_pbc())
        pos = convert(pos, "distance", "ASE", self.units)

        # If necessary, transform the positions to new coordinate system
        if self.coord_transform is not None:
            pos = np.dot(self.coord_transform, pos.transpose())
            pos = pos.transpose()

        # Convert ase_koopmans position matrix to lammps-style position array
        # contiguous in memory
        lmp_positions = list(pos.ravel())

        # Convert that lammps-style array into a C object
        c_double_array = (ctypes.c_double * len(lmp_positions))
        lmp_c_positions = c_double_array(*lmp_positions)
        #        self.lmp.put_coosrds(lmp_c_positions)
        self.lmp.scatter_atoms('x', 1, 3, lmp_c_positions)

    def calculate(self, atoms, properties, system_changes):
        self.propagate(atoms, properties, system_changes, 0)

    def propagate(self, atoms, properties, system_changes, n_steps, dt=None,
                  dt_not_real_time=False, velocity_field=None):
        """"atoms: Atoms object
            Contains positions, unit-cell, ...
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these five: 'positions', 'numbers', 'cell',
            'pbc', 'charges' and 'magmoms'.
        """
        if len(system_changes) == 0:
            return

        self.coord_transform = None

        if not self.started:
            self.start_lammps()
        if not self.initialized:
            self.initialise_lammps(atoms)
        else:  # still need to reset cell
            # Apply only requested boundary condition changes. Note this needs to happen
            # before the call to set_cell since 'change_box' will apply any
            # shrink-wrapping *after* it's updated the cell dimensions
            if 'pbc' in system_changes:
                change_box_str = 'change_box all boundary {}'
                change_box_cmd = change_box_str.format(self.lammpsbc(atoms))
                self.lmp.command(change_box_cmd)

            # Reset positions so that if they are crazy from last
            # propagation, change_box (in set_cell()) won't hang.
            # Could do this only after testing for crazy positions?
            # Could also use scatter_atoms() to set values (requires
            # MPI comm), or extra_atoms() to get pointers to local
            # data structures to zero, but then we would have to be
            # careful with parallelism.
            self.lmp.command("set atom * x 0.0 y 0.0 z 0.0")
            self.set_cell(atoms, change=True)

        if self.parameters.atom_types is None:
            raise NameError("atom_types are mandatory.")

        do_rebuild = (not np.array_equal(atoms.numbers, self.previous_atoms_numbers)
                      or ("numbers" in system_changes))
        if not do_rebuild:
            do_redo_atom_types = not np.array_equal(
                atoms.numbers, self.previous_atoms_numbers)
        else:
            do_redo_atom_types = False

        self.lmp.command('echo none')  # don't echo the atom positions
        if do_rebuild:
            self.rebuild(atoms)
        elif do_redo_atom_types:
            self.redo_atom_types(atoms)
        self.lmp.command('echo log')  # switch back log

        self.set_lammps_pos(atoms)

        if self.parameters.amendments is not None:
            for cmd in self.parameters.amendments:
                self.lmp.command(cmd)

        if n_steps > 0:
            if velocity_field is None:
                vel = convert(atoms.get_velocities(), "velocity", "ASE", self.units)
            else:
                # FIXME: Do we need to worry about converting to lammps units here?
                vel = atoms.arrays[velocity_field]

            # If necessary, transform the velocities to new coordinate system
            if self.coord_transform is not None:
                vel = np.dot(self.coord_transform, vel.T).T

            # Convert ase_koopmans velocities matrix to lammps-style velocities array
            lmp_velocities = list(vel.ravel())

            # Convert that lammps-style array into a C object
            c_double_array = (ctypes.c_double * len(lmp_velocities))
            lmp_c_velocities = c_double_array(*lmp_velocities)
            self.lmp.scatter_atoms('v', 1, 3, lmp_c_velocities)

        # Run for 0 time to calculate
        if dt is not None:
            if dt_not_real_time:
                self.lmp.command('timestep %.30f' % dt)
            else:
                self.lmp.command('timestep %.30f' %
                                 convert(dt, "time", "ASE", self.units))
        self.lmp.command('run %d' % n_steps)

        if n_steps > 0:
            # TODO this must be slower than native copy, but why is it broken?
            pos = np.array(
                [x for x in self.lmp.gather_atoms("x", 1, 3)]).reshape(-1, 3)
            if self.coord_transform is not None:
                pos = np.dot(pos, self.coord_transform)

            # Convert from LAMMPS units to ASE units
            pos = convert(pos, "distance", self.units, "ASE")

            atoms.set_positions(pos)

            vel = np.array(
                [v for v in self.lmp.gather_atoms("v", 1, 3)]).reshape(-1, 3)
            if self.coord_transform is not None:
                vel = np.dot(vel, self.coord_transform)
            if velocity_field is None:
                atoms.set_velocities(convert(vel, 'velocity', self.units,
                                             'ASE'))

        # Extract the forces and energy
        self.results['energy'] = convert(self.lmp.extract_variable('pe', None, 0),
                                         "energy", self.units, "ASE")
        self.results['free_energy'] = self.results['energy']

        stress = np.empty(6)
        stress_vars = ['pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy']

        for i, var in enumerate(stress_vars):
            stress[i] = self.lmp.extract_variable(var, None, 0)

        stress_mat = np.zeros((3, 3))
        stress_mat[0, 0] = stress[0]
        stress_mat[1, 1] = stress[1]
        stress_mat[2, 2] = stress[2]
        stress_mat[1, 2] = stress[3]
        stress_mat[2, 1] = stress[3]
        stress_mat[0, 2] = stress[4]
        stress_mat[2, 0] = stress[4]
        stress_mat[0, 1] = stress[5]
        stress_mat[1, 0] = stress[5]
        if self.coord_transform is not None:
            stress_mat = np.dot(self.coord_transform.T,
                                np.dot(stress_mat, self.coord_transform))
        stress[0] = stress_mat[0, 0]
        stress[1] = stress_mat[1, 1]
        stress[2] = stress_mat[2, 2]
        stress[3] = stress_mat[1, 2]
        stress[4] = stress_mat[0, 2]
        stress[5] = stress_mat[0, 1]

        self.results['stress'] = convert(-stress, "pressure", self.units, "ASE")

        # definitely yields atom-id ordered force array
        f = convert(np.array(self.lmp.gather_atoms("f", 1, 3)).reshape(-1,3),
                    "force", self.units, "ASE")

        if self.coord_transform is not None:
            self.results['forces'] = np.dot(f, self.coord_transform)
        else:
            self.results['forces'] = f.copy()

        # otherwise check_state will always trigger a new calculation
        self.atoms = atoms.copy()

        if not self.parameters.keep_alive:
            self.lmp.close()

    def lammpsbc(self, atoms):
        """
        Determine LAMMPS boundary types base_koopmansd on ASE pbc settings. For non-periodic
        dimensions, if the cell length is finite then fixed BCs ('f') are used; if the
        cell length is approximately zero, shrink-wrapped BCs ('s') are used.
        """
        retval = ''
        pbc = atoms.get_pbc()
        if np.all(pbc):
            retval = 'p p p'
        else:
            cell = atoms.get_cell()
            for i in range(0, 3):
                if pbc[i]:
                    retval += 'p '
                else:
                    # See if we're using indefinite ASE boundaries along this direction
                    if np.linalg.norm(cell[i]) < np.finfo(cell[i][0]).tiny:
                        retval += 's '
                    else:
                        retval += 'f '

        return retval.strip()

    def rebuild(self, atoms):
        try:
            n_diff = len(atoms.numbers) - len(self.previous_atoms_numbers)
        except:
            n_diff = len(atoms.numbers)

        if n_diff > 0:
            if any([("reax/c" in cmd) for cmd in self.parameters.lmpcmds]):
                self.lmp.command("pair_style lj/cut 2.5")
                self.lmp.command("pair_coeff * * 1 1")

                for cmd in self.parameters.lmpcmds:
                    if (("pair_style" in cmd) or ("pair_coeff" in cmd) or
                            ("qeq/reax" in cmd)):
                        self.lmp.command(cmd)

            cmd = "create_atoms 1 random {} 1 NULL".format(n_diff)
            self.lmp.command(cmd)
        elif n_diff < 0:
            cmd = "group delatoms id {}:{}".format(
                len(atoms.numbers) + 1, len(self.previous_atoms_numbers))
            self.lmp.command(cmd)
            cmd = "delete_atoms group delatoms"
            self.lmp.command(cmd)

        self.redo_atom_types(atoms)

    def redo_atom_types(self, atoms):
        current_types = set(
            (i + 1, self.parameters.atom_types[sym]) for i, sym
            in enumerate(atoms.get_chemical_symbols()))

        try:
            previous_types = set(
                (i + 1, self.parameters.atom_types[ase_koopmans_chemical_symbols[Z]])
                for i, Z in enumerate(self.previous_atoms_numbers))
        except:
            previous_types = set()

        for (i, i_type) in current_types - previous_types:
            cmd = "set atom {} type {}".format(i, i_type)
            self.lmp.command(cmd)

        self.previous_atoms_numbers = atoms.numbers.copy()

    def restart_lammps(self, atoms):
        if self.started:
            self.lmp.command("clear")
        # hope there's no other state to be reset
        self.started = False
        self.initialized = False
        self.previous_atoms_numbers = []
        self.start_lammps()
        self.initialise_lammps(atoms)

    def start_lammps(self):
        # Only import lammps when running a calculation
        # so it is not required to use other parts of the
        # module
        from lammps import lammps
        # start lammps process
        if self.parameters.log_file is None:
            cmd_args = ['-echo', 'log', '-log', 'none', '-screen', 'none',
                        '-nocite']
        else:
            cmd_args = ['-echo', 'log', '-log', self.parameters.log_file,
                        '-screen', 'none', '-nocite']

        self.cmd_args = cmd_args

        if self.lmp is None:
            self.lmp = lammps(self.parameters.lammps_name, self.cmd_args,
                              comm=self.parameters.comm)

        # Run header commands to set up lammps (units, etc.)
        for cmd in self.parameters.lammps_header:
            self.lmp.command(cmd)

        for cmd in self.parameters.lammps_header:
            if "units" in cmd:
                self.units = cmd.split()[1]

        if 'lammps_header_extra' in self.parameters:
            if self.parameters.lammps_header_extra is not None:
                for cmd in self.parameters.lammps_header_extra:
                    self.lmp.command(cmd)

        self.started = True

    def initialise_lammps(self, atoms):
        # Initialising commands
        if self.parameters.boundary:
            # if the boundary command is in the supplied commands use that
            # otherwise use atoms pbc
            for cmd in self.parameters.lmpcmds:
                if 'boundary' in cmd:
                    break
            else:
                self.lmp.command('boundary ' + self.lammpsbc(atoms))

        # Initialize cell
        self.set_cell(atoms, change=not self.parameters.create_box)

        if self.parameters.atom_types is None:
            # if None is given, create from atoms object in order of appearance
            s = atoms.get_chemical_symbols()
            _, idx = np.unique(s, return_index=True)
            s_red = np.array(s)[np.sort(idx)].tolist()
            self.parameters.atom_types = {j: i + 1 for i, j in enumerate(s_red)}

        # Initialize box
        if self.parameters.create_box:
            # count number of known types
            n_types = len(self.parameters.atom_types)
            create_box_command = 'create_box {} cell'.format(n_types)
            self.lmp.command(create_box_command)

        # Initialize the atoms with their types
        # positions do not matter here
        if self.parameters.create_atoms:
            self.lmp.command('echo none')  # don't echo the atom positions
            self.rebuild(atoms)
            self.lmp.command('echo log')  # turn back on
        else:
            self.previous_atoms_numbers = atoms.numbers.copy()

        # execute the user commands
        for cmd in self.parameters.lmpcmds:
            self.lmp.command(cmd)

        # Set masses after user commands, e.g. to override
        # EAM-provided masses
        for sym in self.parameters.atom_types:
            if self.parameters.atom_type_masses is None:
                mass = ase_koopmans_atomic_masses[ase_koopmans_atomic_numbers[sym]]
            else:
                mass = self.parameters.atom_type_masses[sym]
            self.lmp.command('mass %d %.30f' % (
                self.parameters.atom_types[sym],
                convert(mass, "mass", "ASE", self.units)))

        # Define force & energy variables for extraction
        self.lmp.command('variable pxx equal pxx')
        self.lmp.command('variable pyy equal pyy')
        self.lmp.command('variable pzz equal pzz')
        self.lmp.command('variable pxy equal pxy')
        self.lmp.command('variable pxz equal pxz')
        self.lmp.command('variable pyz equal pyz')

        # I am not sure why we need this next line but LAMMPS will
        # raise an error if it is not there. Perhaps it is needed to
        # ensure the cell stresses are calculated
        self.lmp.command('thermo_style custom pe pxx emol ecoul')

        self.lmp.command('variable fx atom fx')
        self.lmp.command('variable fy atom fy')
        self.lmp.command('variable fz atom fz')

        # do we need this if we extract from a global ?
        self.lmp.command('variable pe equal pe')

        self.lmp.command("neigh_modify delay 0 every 1 check yes")

        self.initialized = True


# keep this one for the moment being...
def write_lammps_data(filename, atoms, atom_types, comment=None, cutoff=None,
                      molecule_ids=None, charges=None, units='metal'):

    if isinstance(filename, str):
        fh = open(filename, 'w')
    else:
        fh = filename

    if comment is None:
        comment = 'lammpslib autogenerated data file'
    fh.write(comment.strip() + '\n\n')

    fh.write('{0} atoms\n'.format(len(atoms)))
    fh.write('{0} atom types\n'.format(len(atom_types)))

    fh.write('\n')
    cell, coord_transform = convert_cell(atoms.get_cell())
    fh.write('{0:16.8e} {1:16.8e} xlo xhi\n'.format(0.0, cell[0, 0]))
    fh.write('{0:16.8e} {1:16.8e} ylo yhi\n'.format(0.0, cell[1, 1]))
    fh.write('{0:16.8e} {1:16.8e} zlo zhi\n'.format(0.0, cell[2, 2]))
    fh.write('{0:16.8e} {1:16.8e} {2:16.8e} xy xz yz\n'
             ''.format(cell[0, 1], cell[0, 2], cell[1, 2]))

    fh.write('\nMasses\n\n')
    sym_mass = {}
    masses = atoms.get_masses()
    symbols = atoms.get_chemical_symbols()
    for sym in atom_types:
        for i in range(len(atoms)): # TODO: Make this more efficient
            if symbols[i] == sym:
                sym_mass[sym] = convert(masses[i], "mass", "ASE", units)
                break
            else:
                sym_mass[sym] = convert(
                        ase_koopmans_atomic_masses[ase_koopmans_chemical_symbols.index(sym)],
                        "mass", "ASE", units)

    for (sym, typ) in sorted(atom_types.items(), key=operator.itemgetter(1)):
        fh.write('{0} {1}\n'.format(typ, sym_mass[sym]))

    fh.write('\nAtoms # full\n\n')
    if molecule_ids is None:
        molecule_ids = np.zeros(len(atoms), dtype=int)
    if charges is None:
        charges = atoms.get_initial_charges()
    for i, (sym, mol, q, pos) in enumerate(
            zip(symbols, molecule_ids, charges, atoms.get_positions())):
        typ = atom_types[sym]
        fh.write('{0} {1} {2} {3:16.8e} {4:16.8e} {5:16.8e} {6:16.8e}\n'
                 .format(i + 1, mol, typ, q, pos[0], pos[1], pos[2]))

    if isinstance(filename, str):
        fh.close()
