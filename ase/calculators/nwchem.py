"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os
import re
import numpy as np
import time
import collections

from warnings import warn
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem
from ase.calculators.calculator import (FileIOCalculator, Parameters,
                                        all_changes)

# Depending on what type of calculation is being performed, the results can
# have one of several different tag prefixes.
_calc_prefix = ['dft', 'scf', 'mp2', 'ccsd', 'tce', 'pspw', 'band', 'paw']

_etags = [':'.join([prefix, 'energy']) for prefix in _calc_prefix]
_gtags = [':'.join([prefix, 'gradient']) for prefix in _calc_prefix]
_dtags = [':'.join([prefix, 'dipole']) for prefix in _calc_prefix]
_stags = [':'.join([prefix, 'stress']) for prefix in _calc_prefix]

err_mult = 'Multiple {} found in NWchem output file!'

base_keywords = ['geometry']

# Calculator defined keywords that have special treatment
_special_kws = ['center', 'autosym', 'autoz', 'theory', 'basis', 'xc', 'task',
                'pseudopotentials', 'set', 'symmetry']

xc = dict(
        lda='slater pw91lda',
        pbe='xpbe96 cpbe96',
        revpbe='revpbe cpbe96',
        rpbe='rpbe cpbe96',
        pw91='xperdew91 perdew91',
        )

class NWChem(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'stress', 'dipole']
    command = 'nwchem PREFIX.nw > PREFIX.out'

    default_parameters = dict(
        xc='LDA',
        basis='3-21G',
        center=False,
        autosym=False,
        autoz=False,
        )
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='nwchem', atoms=None, command=None, perm=None,
                 scratch=None, task=None, **kwargs):
        """Construct NWchem-calculator object."""
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command, **kwargs)

        if perm is None:
            self.perm = label
        else:
            self.perm = perm
        self.perm = os.path.abspath(self.perm)

        if scratch is None:
            self.scratch = label
        else:
            self.scratch = scratch
        self.scratch = os.path.abspath(self.scratch)

        self.task = task

        p = self.parameters
        if 'theory' in p:
            self.theory = p.theory
        else:
            dft = p.get('dft')
            nwpw = p.get('nwpw')
            if 'xc' in p:
                self.xc = p.get('xc')
                if self.xc in ['ccsd', 'mp2']:
                    self.theory = xc.get(self.xc, self.xc)
                elif self.xc == 'hf':
                    self.theory = 'scf'
                else:
                    if nwpw is not None:
                        self.theory = 'band'
                    else:
                        self.theory = 'dft'
            elif 'tce' in p:
                self.theory = 'tce'
            elif 'ccsd' in p:
                self.theory = 'ccsd'
            elif 'mp2' in p:
                self.theory = 'mp2'
            elif 'scf' in p:
                self.theory = 'hf'
            elif nwpw is not None and 'xc' in nwpw:
                self.theory = 'band'
            else:
                # Do LDA calculation by default in absence of other info
                self.xc = 'lda'
                self.theory = 'dft'

        if self.theory not in ['pspw', 'band', 'paw'] and 'basis' not in p:
            p.basis = '3-21G'

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)

    def _get_geom(self, atoms=None):
        if atoms is None:
            atoms = self.atoms
        outpos = atoms.get_positions()

        geomarg = self.parameters.get('geometry')
        geomline = 'geometry units angstrom'
        if not self.parameters.center:
            geomline += ' nocenter'
        if not self.parameters.autosym:
            geomline += ' noautosym'
        if not self.parameters.autoz:
            geomline += ' noautoz'
        #geomline = 'geometry nocenter noautosym noautoz noprint units angstrom'
        if geomarg is not None:
            geomline += ' ' + geomarg
        geom = [geomline]

        pbc = atoms.pbc

        if np.any(pbc):
            scpos = atoms.get_scaled_positions()
            for i in range(3):
                if pbc[i]:
                    outpos[:, i] = scpos[:, i]

            npbc = pbc.sum()
            cellpars = atoms.get_cell_lengths_and_angles()
            if npbc == 1:
                system = 'polymer'
            elif npbc == 2:
                system = 'surface'
            else:
                system = 'crystal'
            geom.append('  system {} units angstrom'.format(system))
            if npbc == 3:
                geom.append('    lattice_vectors')
                for row in atoms.cell:
                    geom.append('      {:20.16e} {:20.16e} {:20.16e}'.format(*row))
            else:
                if pbc[0]:
                    geom.append('    lat_a {:20.16e}'.format(cellpars[0]))
                if pbc[1]:
                    geom.append('    lat_b {:20.16e}'.format(cellpars[1]))
                if pbc[2]:
                    geom.append('    lat_c {:20.16e}'.format(cellpars[2]))
                if pbc[1] and pbc[2]:
                    geom.append('    alpha {:20.16e}'.format(cellpars[3]))
                if pbc[0] and pbc[2]:
                    geom.append('    beta {:20.16e}'.format(cellpars[4]))
                if pbc[1] and pbc[0]:
                    geom.append('    gamma {:20.16e}'.format(cellpars[5]))
            geom.append('  end')

        for i, atom in enumerate(atoms):
            geom.append('{:>4} {:20.16e} {:20.16e} {:20.16e}'.format(atom.symbol, *outpos[i]))
        symm = self.parameters.get('symmetry')
        if symm is not None:
            geom.append('symmetry {}'.format(symm))
        geom.append('end')

        return geom

    def _get_basis(self):
        if self.theory in ['pspw', 'band', 'paw']:
            return ''
        basis = self.parameters.get('basis')
        if isinstance(basis, str):
            basis = {'*': basis}
        res = ['basis noprint']
        for symbol, ibasis in basis.items():
            res.append('{:>4} library {}'.format(symbol, ibasis))
        res.append('end')

        return res

    def _get_set(self):
        setargs = self.parameters.get('set', dict())
        res = []
        for key, val in setargs.items():
            res.append('set {} {}'.format(key, val))
        return res

    def _format(self, section, params, nindent=0):
        prefix = '  ' * nindent
        if isinstance(params, dict):
            res = [prefix + section]
            for key, val in params.items():
                res += self._format(key, val, nindent+1)
            res.append(prefix + 'end')
            return res
        else:
            return [prefix + '{} {}'.format(section, params)]

    def _get_other(self, params=None, atoms=None):
        if params is None:
            params = self.parameters
        if atoms is None:
            atoms = self.atoms
        res = []

        for kw, block in params.items():
            if kw in _special_kws:
                continue
            #if kw in ['theory', 'geometry', 'basis', 'xc', 'task', 'pseudopotentials', 'set', 'symmetry']:
            #    continue
            res += self._format(kw, block)
        res += self._get_set()
        return res

    def write_input(self, atoms=None, properties=['energy'],
                    system_changes=all_changes):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        os.makedirs(self.perm, exist_ok=True)
        os.makedirs(self.scratch, exist_ok=True)

        if self.task is None:
            if 'forces' in properties or 'stress' in properties:
                task = 'gradient'
            else:
                task = 'energy'
        else:
            task = self.task

        p = self.parameters
        out = ['title "{}"'.format(self.label),
               'permanent_dir {}'.format(self.perm),
               'scratch_dir {}'.format(self.scratch),
               'start {}'.format(self.label),
               '\n'.join(self._get_geom(atoms)),
               '\n'.join(self._get_basis()),
               '\n'.join(self._get_other(atoms=atoms)),
               'task {} {}'.format(self.theory, task),
               'task rtdbprint']

        with open('{}.nw'.format(self.label), 'w') as f:
            f.write('\n\n'.join(out))


    def read_results(self):
        ang_to_Bohr = None
        pos = None
        natoms = None
        energy = None
        forces = None
        stress = None
        dipole = None
        cell = None

        with open(self.label + '.out', 'r') as f:
            lines = f.readlines()

        nlines = len(lines)
        idx = 0
        rtdb_block = False
        while idx < nlines:
            line = lines[idx].strip()
            # Search for the RTDB block at the end of the file
            if not rtdb_block:
                if line.startswith('Contents of RTDB'):
                    rtdb_block = True
                idx += 1
                continue

            # NWChem's builtin Angstrom -> Bohr conversion factor
            if line.startswith('geometry:geometry:angstrom_to_au'):
                ang_to_Bohr = float(lines[idx + 1].split()[0])
                idx += 1

            # System geometry. NWChem lets you keep track of multiple different
            # molecular geometries with different "label"s, the default label
            # being "geometry" (which is why the tag starts geoemtry:geometry:).
            # Currently, we only support this default label.
            elif line.startswith('geometry:geometry:coords'):
                #if pos is not None:
                #    raise RuntimeError(err_mult.format('geometries'))
                ndof = int(line.split()[1][6:].strip('[]'))
                assert ndof % 3 == 0
                if natoms is not None:
                    assert natoms == ndof // 3
                else:
                    natoms = ndof // 3
                pos = []
                nread = 0
                while nread < ndof:
                    idx += 1
                    line = lines[idx].strip().split()
                    nread += len(line)
                    pos += [float(xyz) for xyz in line]
                pos = np.array(pos).reshape((natoms, 3))
                idx += 1

            # "tags" contains the element symbols
            elif line.startswith('geometry:geometry:tags'):
                #if symbols is not None:
                #    raise RuntimeError(err_mult.format('geometry tags'))
                taglen = int(line.split()[1][4:].strip('[]'))
                if natoms is not None:
                    assert natoms == taglen // 2
                else:
                    natoms = taglen // 2
                symbols = []
                idx += 1
                for i in range(natoms):
                    symbols.append(lines[idx + i].strip())
                idx += natoms

            # The unit cell
            elif line.startswith('cell_default:unita_frozen'):
                #if cell is not None:
                #    raise RuntimeError(err_mult.format('unit cells'))
                cell = []
                nread = 0
                while nread < 9:
                    idx += 1
                    line = lines[idx].strip().split()
                    nread += len(line)
                    cell += [float(xyz) for xyz in line]
                cell = np.array(cell).reshape((3, 3))
                idx += 1

            # The energy
            elif any([line.startswith(etag) for etag in _etags]):
                #if energy is not None:
                #    raise RuntimeError(err_mult.format('energies'))
                energy = float(lines[idx + 1].split()[0])
                idx += 2

            # The gradient
            elif any([line.startswith(gtag) for gtag in _gtags]):
                #if forces is not None:
                #    raise RuntimeError(err_mult.format('gradient vectors'))
                ndof = int(line.split()[1][6:].strip('[]'))
                assert ndof % 3 == 0
                if natoms is not None:
                    assert natoms == ndof // 3
                else:
                    natoms = ndof // 3
                grad = []
                nread = 0
                while nread < ndof:
                    idx += 1
                    line = lines[idx].strip().split()
                    nread += len(line)
                    grad += [float(xyz) for xyz in line]
                forces = -np.array(grad).reshape((natoms, 3))
                idx += 1

            # The stress tensor
            elif any([line.startswith(stag) for stag in _stags]):
                #if stress is not None:
                #    raise RuntimeError(err_mult.format('stress tensors'))
                stress = []
                nread = 0
                while nread < 9:
                    idx += 1
                    line = lines[idx].strip().split()
                    nread += len(line)
                    stress += [float(xyz) for xyz in line]
                stress = np.array(stress).reshape((3, 3))
                idx += 1

            # The dipole
            elif any([line.startswith(dtag) for dtag in _dtags]):
                #if dipole is not None:
                #    raise RuntimeError(err_mult.format('dipole moments'))
                dipole = np.array([float(xyz) for xyz in lines[idx + 1].strip().split()])
                idx += 2

            else:
                idx += 1

        # Perform all unit conversions here
        if ang_to_Bohr is None:
            ang_to_Bohr = 1. / Bohr


        if pos is not None:
            self.atoms.positions = pos / ang_to_Bohr
        if cell is not None:
            self.atoms.cell = cell / ang_to_Bohr

        if energy is not None:
            self.results['energy'] = energy * Hartree
        if forces is not None:
            self.results['forces'] = forces * Hartree * ang_to_Bohr
        if stress is not None:
            self.results['stress'] = stress * Hartree * ang_to_Bohr
