"""
This module defines I/O routines with wannier90 files.
For the moment, no attempt has been made to put construct ASE atoms objects correctly.
Instead, everything is stored in ase.calc.parmaeters
"""

import re
import json
import numpy as np
import math
from typing import List
from ase.utils import basestring
from ase.atoms import Atoms
from ase.cell import Cell
from ase.dft.kpoints import BandPath, bandpath
from ase.calculators.wannier90 import Wannier90


def parse_value(value):
    if isinstance(value, list):
        parsed_value = []
        for v in value:
            parsed_value.append(parse_value(v))
    else:
        if isinstance(value, str):
            if value.lower() in ['t', 'true', '.true.']:
                return True
            elif value.lower() in ['f', 'false', '.false.']:
                return False
        try:
            parsed_value = json.loads(value)
        except (TypeError, json.decoder.JSONDecodeError) as e:
            parsed_value = value
    return parsed_value


def write_wannier90_in(fd, atoms):
    """
    Prints out to a given file a wannier90 input file
    """

    settings = {k.lower(): v for k, v in atoms.calc.parameters.items()}

    for kw, opt in settings.items():

        if kw == 'kpoints':
            if np.ndim(opt) == 2:
                rows_str = [' '.join([str(v) for v in row] + [str(1 / len(opt))]) for row in opt]
                fd.write(block_format('kpoints', rows_str))

        elif kw == 'projections':
            rows_str = []
            if 'units' in opt:
                rows_str.append(opt['units'])
            for [site, proj] in opt['sites']:
                if isinstance(site, list):
                    # Interpret as fractional coordinates
                    site_str = 'f=' + ','.join([str(v) for v in site])
                else:
                    # Interpret as atom label
                    site_str = str(site)
                if isinstance(proj, list):
                    proj_str = ';'.join(proj)
                else:
                    proj_str = str(proj)
                rows_str.append(f'{site_str}: {proj_str}')
            fd.write(block_format(kw, rows_str))

        elif kw == 'mp_grid':
            opt_str = ' '.join([str(v) for v in opt])
            fd.write(f'{kw} = {opt_str}\n')

        elif kw == 'kpoint_path':
            if settings.get('bands_plot', False):
                rows_str = []
                for start, end in zip(opt.path[:-1], opt.path[1:]):
                    rows_str.append('')
                    for point in (start, end):
                        rows_str[-1] += ' ' + point + ' '.join([f'{v:9.5f}' for v in opt.special_points[point]])
                fd.write(block_format(kw, rows_str))

        elif isinstance(opt, list):
            if np.ndim(opt) == 1:
                opt = [opt]
            rows_str = [' '.join([str(v) for v in row]) for row in opt]
            fd.write(block_format(kw, rows_str))

        else:
            if isinstance(opt, bool):
                if opt:
                    opt_str = '.true.'
                else:
                    opt_str = '.false.'
            else:
                opt_str = str(opt)
            fd.write(f'{kw} = {opt_str}\n')

    # atoms_frac
    if atoms.has('labels'):
        labels = atoms.get_array('labels')
    else:
        labels = atoms.get_chemical_symbols()
    rows_str = [l + ' ' + ' '.join([str(p) for p in pos]) for l, pos in zip(labels, atoms.get_scaled_positions())]
    fd.write(block_format('atoms_frac', rows_str))

    # unit_cell_cart
    rows_str = ['ang'] + [' '.join([str(v) for v in row]) for row in atoms.cell]
    fd.write(block_format('unit_cell_cart', rows_str))


def block_format(name, rows_str):
    joined_rows = '\n'.join(['  ' + r for r in rows_str])
    return f'\nbegin {name}\n{joined_rows}\nend {name}\n'


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
    kpoint_path = None

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
                    if keyw == 'unit_cell_cart':
                        assert block_lines.pop(0) == ['ang']
                        cell = Cell(parse_value(block_lines))
                    elif keyw == 'atoms_frac':
                        symbols = [l[0] for l in block_lines]
                        scaled_positions = parse_value([l[1:] for l in block_lines])
                    elif keyw == 'projections':
                        block_lines = [[l[0].strip('f=:').split(',')] + l[1:] for l in block_lines]
                        calc.parameters[keyw] = {'sites': parse_value(block_lines)}
                    elif keyw == 'kpoints':
                        calc.parameters[keyw] = np.array([line[:-1] for line in block_lines], dtype=float)
                    elif keyw == 'kpoint_path':
                        kpoint_path = block_lines[0][0]
                        for line, prev_line in zip(block_lines[1:], block_lines[:-1]):
                            if line[0] != prev_line[4]:
                                kpoint_path += prev_line[4] + ','
                            kpoint_path += line[0]
                        kpoint_path += line[4]
                    else:
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
            elif keyw == 'mp_grid':
                calc.parameters[keyw] = [parse_value(v) for v in lsplit[1].split()]
            else:
                calc.parameters[keyw] = parse_value(' '.join(lsplit[1:]))

    # Convert kpoint_path to a BandPath object
    if kpoint_path is not None:
        calc.parameters.kpoint_path = construct_kpoint_path(
            path=kpoint_path, cell=cell, bands_point_num=calc.parameters.get('bands_point_num', 100))

    atoms = Atoms(symbols=symbols, scaled_positions=scaled_positions, cell=cell, calculator=calc)
    atoms.calc.atoms = atoms

    return atoms


def _path_lengths(path: str, cell: Cell, bands_point_num: int) -> List[int]:
    # Construct the metric for reciprocal space (based off Wannier90's "utility_metric" subroutine)
    recip_lat = 2 * math.pi * cell.reciprocal().array
    recip_metric = recip_lat @ recip_lat.T

    # Work out the lengths Wannier90 will assign each path (based off Wannier90's "plot_interpolate_bands" subroutine)
    kpath_pts: List[int] = []
    kpath_len: List[float] = []

    for i, (start, end) in enumerate(zip(path[:-1], path[1:])):
        if start == ',':
            kpath_pts.append(0)
        elif end == ',':
            kpath_pts.append(1)
        else:
            special_points = cell.bandpath().special_points
            vec = special_points[end] - special_points[start]
            kpath_len.append(math.sqrt(vec.T @ recip_metric @ vec))
            if i == 0:
                kpath_pts.append(bands_point_num)
            else:
                kpath_pts.append(int(round(bands_point_num * kpath_len[-1] / kpath_len[0])))
    return kpath_pts


def construct_kpoint_path(path: str, cell: Cell, bands_point_num: int) -> BandPath:
    path_lengths = _path_lengths(path, cell, bands_point_num)

    kpts = []
    for start, end, npoints in zip(path[:-1], path[1:], path_lengths):
        if npoints == 0:
            continue
        bp = bandpath(start + end, cell, npoints + 1).kpts[:-1].tolist()
        kpts += bp
    # Don't forget about the final kpoint
    kpts.append(bp[-1])

    if len(kpts) != sum(path_lengths) + 1:
        raise AssertionError(
            'Did not get the expected number of kpoints; this suggests there is a bug in the code')

    special_points = cell.bandpath().special_points

    return BandPath(cell=cell, kpts=kpts, path=path, special_points=special_points)


def read_wannier90_out(fd):
    """
    Reads wannier90 output files

    Parameters
    ----------
    fd : file|str
        A file like object or filename

    Yields
    ------
    structure : atoms
        An Atoms object with an attached SinglePointCalculator containing
        any parsed results
    """

    if isinstance(fd, basestring):
        fd = open(fd, 'rU')

    flines = fd.readlines()

    structure = Atoms()

    job_done = False
    walltime = None
    convergence = False
    convergence_dis = None
    imre_ratio = []
    centers = []
    spreads = []

    for i, line in enumerate(flines):
        if 'All done' in line:
            job_done = True
        if 'Exiting...' in line and '.nnkp written' in line:
            job_done = True
        if 'Wannierisation convergence criteria satisfied' in line:
            convergence = True
        if 'DISENTANGLE' in line:
            convergence_dis = False
        if 'Disentanglement convergence criteria satisfied' in line:
            convergence_dis = True
        if 'Total Execution Time' in line or 'Time to write kmesh' in line:
            walltime = float(line.split()[-2])
        if 'Maximum Im/Re Ratio' in line:
            imre_ratio.append(float(line.split()[-1]))
        if 'Final State' in line:
            j = 1
            while 'WF centre and spread' in flines[i + j]:
                splitline = [x.rstrip(',') for x in flines[i + j].split()]
                centers.append([float(x) for x in splitline[6:9]])
                spreads.append(float(splitline[-1]))
                j += 1

    calc = Wannier90(atoms=structure)
    calc.results['job done'] = job_done
    calc.results['walltime'] = walltime
    if convergence_dis is not None:
        calc.results['convergence'] = convergence and convergence_dis
    else:
        calc.results['convergence'] = convergence
    calc.results['Im/Re ratio'] = imre_ratio
    calc.results['centers'] = centers
    calc.results['spreads'] = spreads

    structure.calc = calc

    yield structure
