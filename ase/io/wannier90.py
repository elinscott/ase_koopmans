"""
This module defines I/O routines with wannier90 files.
For the moment, no attempt has been made to put construct ASE atoms objects correctly.
Instead, everything is stored in ase.calc.parmaeters
"""

import re
from ase.atoms import Atoms
from ase.calculators.wannier90 import Wannier90
import json

def parse_value(value):
    if isinstance(value, list):
        parsed_value = []
        for v in value:
            parsed_value.append(parse_value(v))
    else:
        if isinstance(value, str):
            if value.lower() in ['t', 'true']:
                return True
            elif value.lower() in ['f', 'false']:
                return False
        try:
           parsed_value = json.loads(value)
        except:
           parsed_value = value
    return parsed_value

def write_wannier90_in(fd, atoms):
    """
    Prints out to a given file a wannier90 input file
    """

    settings = atoms.calc.parameters

    for kw, opt in settings.items():
        if isinstance(opt, list):
            opt_str = '\n'.join([' '.join([str(v) for v in row]) for row in opt])
            fd.write('begin {0}\n{1}\nend {0}\n\n'.format(
                     kw.lower(), opt_str))
        else:
            fd.write('{0} = {1}\n'.format(kw.lower(), opt))


def read_wannier90_in(fd):
    """
    Read a wannier90 input file (the basic format of .cell and .param files)
    and return keyword-value pairs as a dict (values are strings for single
    keywords and lists of strings for blocks).

    Based on ase.io.castep.read_freeform
    """

    filelines = fd.readlines()

    keyw = None
    read_block = False
    block_lines = None

    calc = Wannier90()

    for i, l in enumerate(filelines):

        # Strip all comments, aka anything after a hash
        L = re.split(r'[#!]', l, 1)[0].strip()

        if L == '':
            # Empty line... skip
            continue

        lsplit = re.split(r'\s*[:=]*\s+', L, 1)

        if read_block:
            if lsplit[0].lower() == 'end':
                if len(lsplit) == 1 or lsplit[1].lower() != keyw:
                    raise ValueError(f'Out of place end of block at line {i+1}')
                else:
                    read_block = False
                    calc.parameters[keyw] = parse_value(block_lines)
            else:
                block_lines += [L.split()]
        else:
            # Check the first word

            # Is it a block?
            read_block = (lsplit[0].lower() == 'begin')
            if read_block:
                if len(lsplit) == 1:
                    raise ValueError(f'Unrecognizable block at line {i+1}')
                else:
                    keyw = lsplit[1].lower()
            else:
                keyw = lsplit[0].lower()

            # Now save the value
            if read_block:
                block_lines = []
            else:
                calc.parameters[keyw] = parse_value(' '.join(lsplit[1:]))

    atoms = Atoms(calculator=calc)
    atoms.calc.atoms = atoms

    return atoms
