import os
import re
from subprocess import Popen
from copy import deepcopy

import numpy as np

from ase import Atoms
from ase.utils import workdir
from ase.units import Hartree, Bohr, Debye
from ase.calculators.singlepoint import SinglePointCalculator


def _format_value(val):
    if isinstance(val, bool):
        return '.t.' if val else '.f.'
    return str(val).upper()


def _write_block(name, args):
    out = [' ${}'.format(name.upper())]
    for key, val in args.items():
        out.append('  {}={}'.format(key.upper(), _format_value(val)))
    out.append(' $END')
    return '\n'.join(out)


def _write_geom(atoms):
    out = [' $DATA', atoms.get_chemical_formula(), 'C1']
    for atom in atoms:
        out.append('{} {:>2} {:24.16e} {:24.16e} {:24.16e}'
                   .format(atom.symbol, atom.number, *atom.position))
    out.append(' $END')
    return '\n'.join(out)


def write_gamessus_in(fd, atoms, properties=None, **params):
    params = deepcopy(params)

    if properties is None:
        properties = ['energy']

    contrl = params.pop('contrl', dict())
    if 'runtyp' not in contrl:
        if 'forces' in properties:
            contrl['runtyp'] = 'gradient'
        else:
            contrl['runtyp'] = 'energy'

    out = [_write_block('contrl', contrl)]
    out += [_write_block(*item) for item in params.items()]
    out.append(_write_geom(atoms))
    fd.write('\n\n'.join(out))


_geom_re = re.compile(r'^\s*ATOM\s+ATOMIC\s+COORDINATES')
_atom_re = re.compile(r'^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')
_energy_re = re.compile(r'^\s*FINAL [\S\s]+ ENERGY IS\s+(\S+) AFTER')
_grad_re = re.compile(r'^\s*GRADIENT OF THE ENERGY\s*')
_dipole_re = re.compile(r'^\s+DX\s+DY\s+DZ\s+\/D\/\s+\(DEBYE\)')


def read_gamessus_out(fd):
    atoms = None
    energy = None
    forces = None
    dipole = None
    for line in fd:
        # Geometry
        if _geom_re.match(line):
            fd.readline()
            symbols = []
            pos = []
            while True:
                atom = _atom_re.match(fd.readline())
                if atom is None:
                    break
                symbol, _, x, y, z = atom.groups()
                symbols.append(symbol)
                pos.append(list(map(float, [x, y, z])))
            atoms = Atoms(symbols, np.array(pos) * Bohr)
            continue

        # Energy
        ematch = _energy_re.match(line)
        if ematch is not None:
            energy = float(ematch.group(1)) * Hartree

        # Gradients
        elif _grad_re.match(line):
            for _ in range(3):
                fd.readline()
            grad = []
            while True:
                atom = _atom_re.match(fd.readline())
                if atom is None:
                    break
                grad.append(list(map(float, atom.groups()[2:])))
            forces = -np.array(grad) * Hartree / Bohr
        elif _dipole_re.match(line):
            dipole = np.array(list(map(float, fd.readline().split()[:3])))
            dipole *= Debye

    atoms.calc = SinglePointCalculator(atoms, energy=energy,
                                       forces=forces, dipole=dipole)
    return atoms


def clean_userscr(userscr, prefix):
    for fname in os.listdir(userscr):
        tokens = fname.split('.')
        if tokens[0] == prefix and tokens[-1] != 'bak':
            fold = os.path.join(userscr, fname)
            os.rename(fold, fold + '.bak')


def get_userscr(prefix, command):
    prefix_test = prefix + '_test'
    with workdir(prefix_test, mkdir=True):
        command = command.replace('PREFIX', prefix_test)
        Popen(command, shell=True).wait()

        with open(prefix_test + '.log') as f:
            for line in f:
                if line.startswith('GAMESS supplementary output files'):
                    return ' '.join(line.split(' ')[8:]).strip()

    return None
