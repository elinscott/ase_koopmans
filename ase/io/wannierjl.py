"""
This module defines I/O routines with WannierJL files

In reality, WannierJL does not rely on input and output files, save for a config TOML file for specifying groups
"""

from pathlib import Path

import numpy as np
import toml

from ase import Atoms
from ase.calculators.wannierjl import WannierJL


def write_wannierjl_in(fd, atoms):
    """
    Writes a WannierJL input file (following the toml format)
    """
    # Prepare the dictionary of settings to write to file
    dct = {}

    def sanitize(obj):
        if isinstance(obj, Path):
            obj = str(obj)
        elif isinstance(obj, list):
            obj = [sanitize(x) for x in obj]
        return obj
    dct['groups'] = {k: sanitize(v) for k, v in atoms.calc.parameters.items()
                     if k in ['indices', 'outdirs']}

    # Write the wjl input file
    toml.dump(dct, fd)
    return


def read_wannierjl_in(fd):
    # Create the calculator
    calc = WannierJL()

    # Create the Atoms object and link the calculator
    atoms = Atoms(calculator=calc)
    calc.atoms = atoms

    # Load the calculator parameters from file
    dct = toml.load(fd)

    calc.parameters = dct['groups']

    # Return the atoms object
    return atoms


def read_wannierjl_out(fd):
    """
    Dummy function that pretends to read WannierJL output files

    Parameters
    ----------
    fd : file|str
        A file like object or filename

    Yields
    ------
    structure : atoms
        An Atoms object with an attached WannierJL Calculator containing
        any parsed results
    """

    structure = Atoms()
    calc = WannierJL(atoms=structure)
    structure.calc = calc

    flines = fd.readlines()

    job_done = False
    ngroups = 1
    wannier_tables = []
    for i, line in enumerate(flines):
        if 'Model will be split into' in line:
            ngroups = int(line.split()[-2])
        if line.strip().startswith('WF     center'):
            i_end = flines[i:].index(
                [l for l in flines[i:] if 'Sum spread' in l][0])
            wannier_tables.append(
                [[float(x) for x in l.strip().split()[1:]] for l in flines[i + 1:i + i_end]])

    # Job is done if we have ngroups + 1 tables
    if len(wannier_tables) > ngroups:
        job_done = True
        # Only read the last ngroups tables, concatenated together
        wannier_tables = np.concatenate([np.array(x) for x in wannier_tables[-ngroups:]])
        calc.results['centers'] = wannier_tables[:, :-1]
        calc.results['spreads'] = wannier_tables[:, -1]
    calc.results['job_done'] = job_done

    yield structure
