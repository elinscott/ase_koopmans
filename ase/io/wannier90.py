"""
This module defines I/O routines with wannier90 files.
For the moment, no attempt has been made to put construct ASE atoms objects correctly.
Instead, everything is stored in ase.calc.parmaeters
"""

import re
import json
import numpy as np
import math
import warnings
from typing import List, Dict, Any, Union
from ase.utils import basestring
from ase.atoms import Atoms
from ase.cell import Cell
from ase.geometry import wrap_positions
from ase.dft.kpoints import BandPath, bandpath
from ase.calculators.wannier90 import Wannier90


def list_to_formatted_str(values: List[int]) -> str:
    # Converts a list of integers into the format expected by Wannier90
    # e.g. list_to_formatted_str([1, 2, 3, 4, 5, 7]) = "1-5,7"
    if len(values) == 0:
        raise ValueError('list_to_formatted_str() should not be given an empty list')
    assert all(a > b for a, b in zip(values[1:], values[:-1])), 'values must be monotonically increasing'
    indices: List[Union[int, None]] = [None]
    indices += [i + 1 for i in range(len(values) - 1) if values[i + 1] != values[i] + 1]
    indices += [None]
    sectors = [values[slice(a, b)] for a, b in zip(indices[:-1], indices[1:])]
    out = []
    for sector in sectors:
        if len(sector) == 1:
            out.append(str(sector[0]))
        else:
            out.append(f'{sector[0]}-{sector[-1]}')
    return ','.join(out)


def formatted_str_to_list(string: str) -> List[int]:
    # Performs the inverse of list_to_formatted_str
    out = []
    for section in string.split(','):
        if '-' in section:
            out += list(range(int(section.split('-')[0]), int(section.split('-')[1]) + 1))
        else:
            out.append(int(section))
    return out


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
            rows_str = [proj_dict_to_string(dct, atoms) for dct in opt]
            fd.write(block_format(kw, rows_str))

        elif kw == 'mp_grid':
            opt_str = ' '.join([str(v) for v in opt])
            fd.write(f'{kw} = {opt_str}\n')

        elif kw == 'kpoint_path':
            if settings.get('bands_plot', False):
                rows_str = []
                path_list = path_str_to_list(opt.path, opt.special_points)
                for start, end in zip(path_list[:-1], path_list[1:]):
                    if ',' in [start, end]:
                        continue
                    rows_str.append('')
                    for point in (start, end):
                        rows_str[-1] += ' ' + point + ' '.join([f'{v:9.5f}' for v in opt.special_points[point]])
                fd.write(block_format(kw, rows_str))

        elif isinstance(opt, list) and kw != 'exclude_bands':
            if np.ndim(opt) == 1:
                opt = [opt]
            rows_str = [' '.join([str(v) for v in row]) for row in opt]
            fd.write(block_format(kw, rows_str))

        else:
            if kw == 'exclude_bands' and isinstance(opt, list):
                opt_str = list_to_formatted_str(opt)
            elif isinstance(opt, bool):
                if opt:
                    opt_str = '.true.'
                else:
                    opt_str = '.false.'
            else:
                opt_str = str(opt)
            fd.write(f'{kw} = {opt_str}\n')

    # atoms_frac
    if atoms.has('tags'):
        labels = [s + str(t) if t != 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
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
                        symbols = []
                        tags = []
                        for l in block_lines:
                            symbols.append(''.join([c for c in l[0] if c.isalpha()]))
                            tag = ''.join([c for c in l[0] if c.isnumeric()])
                            if tag:
                                tags.append(int(tag))
                            else:
                                tags.append(0)
                        scaled_positions = parse_value([l[1:] for l in block_lines])
                    elif keyw == 'projections':
                        calc.parameters[keyw] = [proj_string_to_dict(' '.join(line)) for line in block_lines]
                    elif keyw == 'kpoints':
                        calc.parameters[keyw] = np.array([line[:-1] for line in block_lines], dtype=float)
                    elif keyw == 'kpoint_path':
                        if len(block_lines) == 1:
                            kpoint_path = block_lines[0][0] + block_lines[0][4]
                        else:
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
            elif keyw == 'exclude_bands':
                calc.parameters[keyw] = formatted_str_to_list(' '.join(lsplit[1:]))
            else:
                calc.parameters[keyw] = parse_value(' '.join(lsplit[1:]))

    # Convert kpoint_path to a BandPath object
    if kpoint_path is not None:
        calc.parameters.kpoint_path = construct_kpoint_path(
            path=kpoint_path, cell=cell, bands_point_num=calc.parameters.get('bands_point_num', 100))

    atoms = Atoms(symbols=symbols, scaled_positions=scaled_positions, cell=cell, calculator=calc, tags=tags)
    atoms.calc.atoms = atoms

    return atoms


def proj_dict_to_string(dct, atoms):
    site = []
    site_outside_pc = False
    if 'csite' in dct:
        coords = np.array(dct['csite'])
        # Check coords are inside the cell
        [wrapped_coords] = wrap_positions([coords], atoms.cell)
        if np.linalg.norm(wrapped_coords - coords) > 1e-3:
            site_outside_pc = True
            coords = wrapped_coords
        site.append('c=' + ','.join([str(c) for c in coords]))
    elif 'fsite' in dct:
        coords = dct['fsite']
        # Check coords are inside the cell
        if any([x // 1 != 0 for x in coords]):
            site_outside_pc = True
            coords = [x % 1 for x in coords]
        site.append('f=' + ','.join([str(c) for c in coords]))
    elif 'site' in dct:
        site.append(dct['site'])
    else:
        raise ValueError('w90 projections block is missing a "site" entry')

    if site_outside_pc:
        warnings.warn('A projection site lies outside the primitive cell. It has been wrapped.')

    if 'ang_mtm' not in dct:
        raise ValueError('w90 projections block is missing an "ang_mtm" entry')
    site.append(dct['ang_mtm'])

    if 'zaxis' in dct:
        site.append('z=' + ','.join(str(x) for x in dct['zaxis']))

    if 'xaxis' in dct:
        site.append('x=' + ','.join(str(x) for x in dct['zaxis']))

    if 'radial' in dct:
        site.append(f'r={dct["radial"]}')

    if 'zona' in dct:
        site.append(f'zona={dct["zona"]}')

    return ':'.join(site)


def proj_string_to_dict(string):
    # Converts a w90-formatted projector string into a dictionary
    proj = {}
    [site, ang_mtm] = string.split(':')[:2]

    # Storing site
    if '=' in site:
        site_vec = parse_value(site.split('=')[1].split(','))
        if site.startswith('c'):
            proj['csite'] = site_vec
        else:
            proj['fsite'] = site_vec
    else:
        proj['site'] = site

    # Storing ang_mtm
    proj['ang_mtm'] = ang_mtm

    # Storing other optional arguments
    for setting in string.split(':')[2:]:
        k, v = [x.strip() for x in setting.split('=')]

        k_to_label = {'z': 'zaxis', 'x': 'xaxis', 'r': 'radial', 'zona': 'zona'}
        assert k in k_to_label, f'Failed to parse {k} in the projection {string}'

        if ',' in v:
            v = v.split(',')

        proj[k_to_label[k]] = parse_value(v)

    return proj


def _path_lengths(path: str, cell: Cell, bands_point_num: int) -> List[int]:
    # Construct the metric for reciprocal space (based off Wannier90's "utility_metric" subroutine)
    recip_lat = 2 * math.pi * cell.reciprocal().array
    recip_metric = recip_lat @ recip_lat.T

    # Work out the lengths Wannier90 will assign each path (based off Wannier90's "plot_interpolate_bands" subroutine)
    kpath_pts: List[int] = []
    kpath_len: List[float] = []
    special_points = cell.bandpath().special_points

    path_list = path_str_to_list(path, special_points)

    for i, (start, end) in enumerate(zip(path_list[:-1], path_list[1:])):
        if start == ',':
            kpath_pts.append(0)
        elif end == ',':
            kpath_pts.append(1)
        else:
            vec = special_points[end] - special_points[start]
            kpath_len.append(math.sqrt(vec.T @ recip_metric @ vec))
            if i == 0:
                kpath_pts.append(bands_point_num)
            else:
                kpath_pts.append(int(round(bands_point_num * kpath_len[-1] / kpath_len[0])))
    return kpath_pts


def path_str_to_list(path: str, special_points: Dict[str, str]) -> List[str]:
    # Breaks down a k-point path specified by a string into a list of special points e.g.
    # 'GLB1,BZGX,QFP1Z,LP' -> ['G', 'L', 'B1', ',', 'B', 'Z', 'G', 'X', ',', 'Q', 'F', 'P1', 'Z', ',', 'L', 'P']
    path_list = []
    while len(path) > 0:
        found = False
        for option in list(special_points.keys()) + [',']:
            if path.endswith(option):
                path = path[:-len(option)]
                path_list.append(option)
                found = True
                break
        if not found:
            raise ValueError(f'Could not deconstruct {path} into individual high-symmetry points')
    return path_list[::-1]


def construct_kpoint_path(path: str, cell: Cell, bands_point_num: int) -> BandPath:
    path_lengths = _path_lengths(path, cell, bands_point_num)
    special_points = cell.bandpath().special_points

    path_list = path_str_to_list(path, special_points)

    kpts = []
    for start, end, npoints in zip(path_list[:-1], path_list[1:], path_lengths):
        if start == ',':
            pass
        elif end == ',':
            kpts.append(special_points[start].tolist())
        else:
            bp = bandpath(start + end, cell, npoints + 1)
            kpts += bp.kpts[:-1].tolist()
    # Don't forget about the final kpoint
    kpts.append(bp.kpts[-1].tolist())

    if len(kpts) != sum(path_lengths) + 1:
        raise AssertionError(
            'Did not get the expected number of kpoints; this suggests there is a bug in the code')

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
                splitline = [x.rstrip(',') for x in re.sub('[(),]', ' ', flines[i + j]).split()]
                centers.append([float(x) for x in splitline[5:8]])
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


def num_wann_lookup(proj):
    database = {'l=0': 1, 's': 1,
                'l=1': 3, 'p': 3, 'pz': 1, 'px': 1, 'py': 1,
                'l=2': 5, 'd': 5, 'dxy': 1, 'dxz': 1, 'dyz': 1, 'dx2-y2': 1, 'dz2': 1,
                'l=3': 7, 'f': 7, 'fz3': 1, 'fxz2': 1, 'fyz2': 1, 'fz(x2-y2)': 1, 'fxyz': 1, 'fx(x2-3y2)': 1,
                'fy(3x2-y2)': 1,
                'l=-1': 2, 'sp': 2, 'sp-1': 1, 'sp-2': 1,
                'l=-2': 3, 'sp2': 3, 'sp2-1': 1, 'sp2-2': 1, 'sp2-3': 1,
                'l=-3': 4, 'sp3': 4, 'sp3-1': 1, 'sp3-2': 1, 'sp3-3': 1, 'sp3-4': 1,
                'l=-4': 5, 'sp3d': 5, 'sp3d-1': 1, 'sp3d-2': 1, 'sp3d-3': 1, 'sp3d-4': 1, 'sp3d-5': 1,
                'l=-5': 6, 'sp3d2': 6, 'sp3d2-1': 1, 'sp3d2-2': 1, 'sp3d2-3': 1, 'sp3d2-4': 1, 'sp3d2-5': 1,
                'sp3d2-6': 1}
    if proj in database:
        return database[proj]
    elif proj.split(',')[0] in database:
        l, mr = proj.split(',')
        if l in database and '=' in mr:
            _, val = mr.split('=')
            return len(formatted_str_to_list(val))

    raise NotImplementedError(f'Unrecognised Wannier90 projection {proj}')


def num_wann_from_projections(projections: List[Dict[str, Any]], atoms: Atoms):
    # Works out the value of 'num_wann' based on the 'projections' block
    num_wann = 0
    for proj in projections:
        if 'site' in proj:
            if len(set(atoms.get_tags())) > 0:
                labels = [s + str(t) if t != 0 else s for s, t in zip(atoms.symbols, atoms.get_tags())]
            else:
                labels = [s for s in atoms.symbols]
            num_sites = labels.count(proj['site'])
        else:
            num_sites = 1
        ang_mtms = [p.replace(' ', '') for p in proj['ang_mtm'].split(';')]
        num_wann += num_sites * sum([num_wann_lookup(ang_mtm) for ang_mtm in ang_mtms])
    return num_wann


