from __future__ import print_function
# lammps.py (2011/03/29)
# An ASE calculator for the LAMMPS classical MD code available from
#       http://lammps.sandia.gov/
# The environment variable LAMMPS_COMMAND must be defined to point to the
# LAMMPS binary.
#
# Copyright (C) 2009 - 2011 Joerg Meyer, joerg.meyer@ch.tum.de
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA or see <http://www.gnu.org/licenses/>.


import os
import shutil
import shlex
from subprocess import Popen, PIPE
from threading import Thread
from re import compile as re_compile, IGNORECASE
from tempfile import mkdtemp, NamedTemporaryFile, mktemp as uns_mktemp
import numpy as np
from ase import Atoms
from ase.parallel import paropen
from ase.calculators.calculator import Calculator
from ase.calculators.calculator import all_changes
from ase.utils import basestring
from ase.io.lammpsdata import write_lammps_data
from ase.io.lammpsrun import read_lammps_dump
from ase.calculators.lammps import Prism, write_lammps_in, CALCULATION_END_MARK, convert

__all__ = ['LAMMPS']


class LAMMPS(Calculator):
    name = 'lammpsrun'

    implemented_properties = ['energy', 'forces', 'stress']

    default_parameters = dict(
        units='metal',
        atom_style='atomic',
        dump_period=1,
        specorder=None,
        always_triclinic=False,
        verbose=False,
        thermo_args=['step', 'temp', 'press', 'cpu',
                     'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                     'ke', 'pe', 'etotal',
                     'vol', 'lx', 'ly', 'lz', 'atoms']
        )

    # legacy parameter persist, when the 'parameters' is manipulated from the
    # outside.  All others are rested to the default value
    legacy_parameters=['specorder',
                       'dump_period',
                       '_custom_thermo_args',
                       ]
    
    def __init__(self, label='lammps', tmp_dir=None, 
                 parameters=default_parameters, files=[], always_triclinic=False,
                 keep_alive=True, keep_tmp_files=False,
                 no_data_file=False,
                 **kwargs):
        """The LAMMPS calculators object

        files: list
            Short explanation XXX
        parameters: dict
            Short explanation XXX
        keep_tmp_files: bool
            Retain any temporary files created. Mostly useful for debugging.
        tmp_dir: str
            path/dirname (default None -> create automatically).
            Explicitly control where the calculator object should create
            its files. Using this option implies 'keep_tmp_files'
        no_data_file: bool
            Controls whether an explicit data file will be used for feeding
            atom coordinates into lammps. Enable it to lessen the pressure on
            the (tmp) file system. THIS OPTION MIGHT BE UNRELIABLE FOR CERTAIN
            CORNER CASES (however, if it fails, you will notice...).
        keep_alive: bool
            When using LAMMPS as a spawned subprocess, keep the subprocess
            alive (but idling when unused) along with the calculator object.
        always_triclinic: bool
            Force use of a triclinic cell in LAMMPS, even if the cell is
            a perfect parallelepiped.
        """

        Calculator.__init__(self, label=label, **kwargs)

        self.parameters.update(parameters)
        self.prism = None
        self.files = files
        self.calls = 0
        self.forces = None
        # thermo_content contains data "written by" thermo_style.
        # It is a list of dictionaries, each dict (one for each line
        # printed by thermo_style) contains a mapping between each
        # custom_thermo_args-argument and the corresponding
        # value as printed by lammps. thermo_content will be
        # re-populated by the read_log method.
        self.thermo_content = []

        
        self.keep_alive = keep_alive
        self.no_data_file = no_data_file
        self.keep_tmp_files = keep_tmp_files

        # if True writes velocities from atoms.get_velocities() to LAMMPS input
        self.write_velocities = False

        # file object, if is not None the trajectory will be saved in it
        self.trajectory_out = None

        # period of system snapshot saving (in MD steps)
        parameters['dump_period'] = 1
        
        if not hasattr(self.parameters, 'always_triclinic'):
            self.parameters['always_triclinic'] = always_triclinic
        if not hasattr(self.parameters, 'keep_tmp_files'):
            self.parameters['verbose'] = keep_tmp_files

        if tmp_dir is not None:
            # If tmp_dir is pointing somewhere, don't remove stuff!
            self.keep_tmp_files = True
        self._lmp_handle = None        # To handle the lmp process

        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='LAMMPS-')
        else:
            self.tmp_dir = os.path.realpath(tmp_dir)
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir, 0o755)

        for f in files:
            shutil.copy(f, os.path.join(self.tmp_dir, os.path.basename(f)))

    def __setattr__(self, key, value):
        """Old LAMMPSRUN allows it to just override the parameters
        dictionary. "Modern" ase calculators can assume that default
        parameters are always set, overrides of the
        'parameters'-dictionary have to be caught and the default
        parameters need to be added first.
        """
        # !TODO: remove and break somebody's code (e.g. the test example)
        if key == 'parameters' and value is not None:
            temp_dict = self.get_default_parameters()
            if self.parameters:
                for l_key in self.legacy_parameters:
                    try:
                        temp_dict[l_key] = self.parameters[l_key]
                    except KeyError:
                        pass
            temp_dict.update(value)
            value = temp_dict
        Calculator.__setattr__(self, key, value)

    def clean(self, force=False):

        self._lmp_end()

        if not self.parameters['keep_tmp_files']:
            shutil.rmtree(self.tmp_dir)

    def check_state(self, atoms, tol=1.0e-4):
        # differenct convention for unit-cell and limit precision in
        # LAMMPS-input file will lead to small rounding errors
        return Calculator.check_state(self, atoms, tol)

    def calculate(self, atoms=None,
                  properties=['energy', 'forces', 'stress', 'energies'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        pbc = self.atoms.get_pbc()
        if all(pbc):
            cell = self.atoms.get_cell()
        elif not any(pbc):
            # large enough cell for non-periodic calculation -
            # LAMMPS shrink-wraps automatically via input command
            #       "periodic s s s"
            # below
            cell = 2 * np.max(np.abs(self.atoms.get_positions())) * np.eye(3)
        else:
            print("WARNING: semi-periodic ASE cell detected - translation")
            print("         to proper LAMMPS input cell might fail")
            cell = self.atoms.get_cell()
        self.prism = Prism(cell)
        self.run()

        tc = self.thermo_content[-1]
        
        self.results['energy'] = convert(tc['pe'],'energy',
                                         self.parameters['units'], 'ASE')
        self.results['forces'] = self.forces.copy()
        stress = np.array([tc[i] for i in ('pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy')])
        self.results['stress'] = convert(stress, 'pressure',
                                         self.parameters['units'], 'ASE')

    def _lmp_alive(self):
        # Return True if this calculator is currently handling a running
        # lammps process
        return self._lmp_handle and not isinstance(
            self._lmp_handle.poll(), int)

    def _lmp_end(self):
        # Close lammps input and wait for lammps to end. Return process
        # return value
        if self._lmp_alive():
            self._lmp_handle.stdin.close()
            return self._lmp_handle.wait()

    def run(self, set_atoms=False):
        """Method which explicitly runs LAMMPS."""

        self.calls += 1

        # set LAMMPS command from environment variable
        if 'LAMMPS_COMMAND' in os.environ:
            lammps_cmd_line = shlex.split(os.environ['LAMMPS_COMMAND'],
                                          posix=(os.name == 'posix'))
            if len(lammps_cmd_line) == 0:
                self.clean()
                raise RuntimeError('The LAMMPS_COMMAND environment variable '
                                   'must not be empty')
            # want always an absolute path to LAMMPS binary when calling from
            # self.dir
            lammps_cmd_line[0] = os.path.abspath(lammps_cmd_line[0])

        else:
            self.clean()
            raise RuntimeError(
                'Please set LAMMPS_COMMAND environment variable')
        if 'LAMMPS_OPTIONS' in os.environ:
            lammps_options = shlex.split(os.environ['LAMMPS_OPTIONS'],
                                         posix=(os.name == 'posix'))
        else:
            lammps_options = shlex.split('-echo log -screen none',
                                         posix=(os.name == 'posix'))

        # change into subdirectory for LAMMPS calculations
        cwd = os.getcwd()
        os.chdir(self.tmp_dir)

        # setup file names for LAMMPS calculation
        label = '{0}{1:>06}'.format(self.label, self.calls)
        lammps_in = uns_mktemp(prefix='in_' + label, dir=self.tmp_dir)
        lammps_log = uns_mktemp(prefix='log_' + label, dir=self.tmp_dir)
        lammps_trj_fd = NamedTemporaryFile(
            prefix='trj_' + label, suffix='.bin', dir=self.tmp_dir,
            delete=(not self.keep_tmp_files))
        lammps_trj = lammps_trj_fd.name
        if self.no_data_file:
            lammps_data = None
        else:
            lammps_data_fd = NamedTemporaryFile(
                prefix='data_' + label, dir=self.tmp_dir,
                delete=(not self.keep_tmp_files))
            write_lammps_data(lammps_data_fd, self.atoms,
                              specorder=self.parameters['specorder'],
                              force_skew=self.parameters['always_triclinic'],
                              prismobj=self.prism)
            lammps_data = lammps_data_fd.name
            lammps_data_fd.flush()

        # see to it that LAMMPS is started
        if not self._lmp_alive():
            # Attempt to (re)start lammps
            self._lmp_handle = Popen(
                lammps_cmd_line + lammps_options + ['-log', '/dev/stdout'],
                stdin=PIPE, stdout=PIPE)
        lmp_handle = self._lmp_handle

        # Create thread reading lammps stdout (for reference, if requested,
        # also create lammps_log, although it is never used)
        if self.keep_tmp_files:
            lammps_log_fd = open(lammps_log, 'wb')
            fd = SpecialTee(lmp_handle.stdout, lammps_log_fd)
        else:
            fd = lmp_handle.stdout
        thr_read_log = Thread(target=self.read_lammps_log, args=(fd,))
        thr_read_log.start()

        # write LAMMPS input (for reference, also create the file lammps_in,
        # although it is never used)
        if self.keep_tmp_files:
            lammps_in_fd = open(lammps_in, 'wb')
            fd = SpecialTee(lmp_handle.stdin, lammps_in_fd)
        else:
            fd = lmp_handle.stdin
        write_lammps_in(lammps_in=fd,
                        parameters=self.parameters,
                        atoms=self.atoms,
                        prismobj=self.prism,
                        lammps_trj=lammps_trj,
                        lammps_data=lammps_data)

        if self.keep_tmp_files:
            lammps_in_fd.close()

        # Wait for log output to be read (i.e., for LAMMPS to finish)
        # and close the log file if there is one
        thr_read_log.join()
        if self.keep_tmp_files:
            lammps_log_fd.close()

        if not self.keep_alive:
            self._lmp_end()

        exitcode = lmp_handle.poll()
        if exitcode and exitcode != 0:
            cwd = os.getcwd()
            raise RuntimeError('LAMMPS exited in {} with exit code: {}.'
                               ''.format(cwd, exitcode))

        # A few sanity checks
        if len(self.thermo_content) == 0:
            raise RuntimeError('Failed to retrieve any thermo_style-output')
        if int(self.thermo_content[-1]['atoms']) != len(self.atoms):
            # This obviously shouldn't happen, but if prism.fold_...() fails,
            # it could
            raise RuntimeError('Atoms have gone missing')

        trj_atoms = read_lammps_dump(infileobj=lammps_trj,
                                     order=False,
                                     index=-1,
                                     prismobj=self.prism,
                                     specorder=self.parameters['specorder'])

        if set_atoms:
            self.atoms = trj_atoms.copy()

        self.forces = trj_atoms.get_forces()
        # !TODO: trj_atoms is only the last snapshot of the system; Is it
        #        desireable to save also the inbetween steps?
        if self.trajectory_out is not None:
            # !TODO: is it advisable to create here temporary atoms-objects
            self.trajectory_out.write(trj_atoms)

        lammps_trj_fd.close()
        if not self.no_data_file:
            lammps_data_fd.close()

        os.chdir(cwd)

    def read_lammps_log(self, lammps_log=None, PotEng_first=False):
        """Method which reads a LAMMPS output log file."""

        if lammps_log is None:
            lammps_log = self.label + '.log'

        if isinstance(lammps_log, basestring):
            f = paropen(lammps_log, 'wb')
            close_log_file = True
        else:
            # Expect lammps_in to be a file-like object
            f = lammps_log
            close_log_file = False
            
        # read_log depends on that the first (three) thermo_style custom args
        # can be capitilized and matched against the log output. I.e.
        # don't use e.g. 'ke' or 'cpu' which are labeled KinEng and CPU.
        _custom_thermo_mark = ' '.join([x.capitalize() for x in
                                        self.parameters['thermo_args'][0:3]])

        # !TODO: regex-magic necessary?
        # Match something which can be converted to a float
        f_re = r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?|nan|inf))'
        n_args = len(self.parameters['thermo_args'])
        # Create a re matching exactly N white space separated floatish things
        _custom_thermo_re = re_compile(
            r'^\s*' + r'\s+'.join([f_re] * n_args) + r'\s*$', flags=IGNORECASE)
            
        thermo_content = []
        line = f.readline().decode('utf-8')
        while line and line.strip() != CALCULATION_END_MARK:
            # get thermo output
            if line.startswith(_custom_thermo_mark):
                m = True
                while m:
                    line = f.readline().decode('utf-8')
                    m = _custom_thermo_re.match(line)
                    if m:
                        # create a dictionary between each of the
                        # thermo_style args and it's corresponding value
                        thermo_content.append(
                            dict(zip(self.parameters['thermo_args'],
                                     map(float, m.groups()))))
            else:
                line = f.readline().decode('utf-8')

        if close_log_file:
            f.close()

        self.thermo_content = thermo_content


class SpecialTee(object):
    """A special purpose, with limited applicability, tee-like thing.

    A subset of stuff read from, or written to, orig_fd,
    is also written to out_fd.
    It is used by the lammps calculator for creating file-logs of stuff
    read from, or written to, stdin and stdout, respectively.
    """

    def __init__(self, orig_fd, out_fd):
        self._orig_fd = orig_fd
        self._out_fd = out_fd
        self.name = orig_fd.name

    def write(self, data):
        self._orig_fd.write(data)
        self._out_fd.write(data)
        self.flush()

    def read(self, *args, **kwargs):
        data = self._orig_fd.read(*args, **kwargs)
        self._out_fd.write(data)
        return data

    def readline(self, *args, **kwargs):
        data = self._orig_fd.readline(*args, **kwargs)
        self._out_fd.write(data)
        return data

    def readlines(self, *args, **kwargs):
        data = self._orig_fd.readlines(*args, **kwargs)
        self._out_fd.write(''.join(data))
        return data

    def flush(self):
        self._orig_fd.flush()
        self._out_fd.flush()


if __name__ == '__main__':
    pair_style = 'eam'
    Pd_eam_file = 'Pd_u3.eam'
    pair_coeff = ['* * ' + Pd_eam_file]
    parameters = {'pair_style': pair_style, 'pair_coeff': pair_coeff}
    my_files = [Pd_eam_file]
    calc = LAMMPS(parameters=parameters, files=my_files)
    a0 = 3.93
    b0 = a0 / 2.0
    if True:
        bulk = Atoms(
            ['Pd'] * 4,
            positions=[(0, 0, 0), (b0, b0, 0), (b0, 0, b0), (0, b0, b0)],
            cell=[a0] * 3, pbc=True)
        # test get_forces
        print('forces for a = {0}'.format(a0))
        print(calc.get_forces(bulk))
        # single points for various lattice constants
        bulk.set_calculator(calc)
        for i in range(-5, 5, 1):
            a = a0 * (1 + i / 100.0)
            bulk.set_cell([a] * 3)
            print('a : {0} , total energy : {1}'.format(
                a, bulk.get_potential_energy()))

    calc.clean()
