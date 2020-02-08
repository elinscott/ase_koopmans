"""GPUMD input file parser.

Implemented:
* Input file (xyz.in)

"""

import numpy as np

from ase import io
from ase.build import bulk
from ase.io.gpumd import load_xyz_input_gpumd


# This file is parsed correctly by GPUMD, since it include
# among the examples distributed with the package, i.e.
# GPUMD/examples/ex5/xyz.in
gpumd_input_text = """
15 1 0.1 0 0 1
1 1 1 16.2 16.2 16.2
0 0.0993428 -0.0276529 0.129538 26.9815 0
0 0.304606 1.97817 1.97817 26.9815 1
0 2.34084 0.153487 1.93111 26.9815 2
0 2.13351 1.93232 -0.093146 26.9815 3
0 0.0483925 -0.382656 3.70502 26.9815 4
0 -0.112458 1.82243 6.13785 26.9815 5
0 1.8434 -0.282461 6.36813 26.9815 6
0 1.97984 2.03851 3.76505 26.9815 7
0 -0.108877 0.0221845 7.8698 26.9815 8
0 0.0751396 1.90487 10.0667 26.9815 9
0 1.90466 0.370456 10.1223 26.9815 10
0 1.81346 2.18951 7.85583 26.9815 11
0 0.0417727 -0.391934 11.8844 26.9815 12
0 0.0393722 2.17269 14.2093 26.9815 13
0 2.00187 -0.0602207 13.8793 26.9815 14
"""


def test_read_gpumd_input():
    """Read GPUMD input file."""
    with open('xyz.in', 'w') as gpumd_input_f:
        gpumd_input_f.write(gpumd_input_text)

    # Test when specifying the species-type map
    species_types = {'Al': 0}
    gpumd_input_atoms = io.read('xyz.in', format='gpumd',
                                species_types=species_types)
    assert len(gpumd_input_atoms) == 15
    assert all(s == 'Al' for s in gpumd_input_atoms.get_chemical_symbols())
    assert all(gpumd_input_atoms.get_pbc())
    assert len(gpumd_input_atoms.info) == len(gpumd_input_atoms)
    assert all(np.array_equal(
        gpumd_input_atoms.info[i]['groups'], np.array([i])) for i in
        range(len(gpumd_input_atoms)))

    # Test without specifying the species-type map
    gpumd_input_atoms = io.read('xyz.in', format='gpumd')
    assert all(s == 'Al' for s in gpumd_input_atoms.get_chemical_symbols())

    # Test when specifying the isotope masses
    isotope_masses = {'Al': [26.9815]}
    gpumd_input_atoms = io.read('xyz.in', format='gpumd',
                                isotope_masses=isotope_masses)
    assert all(s == 'Al' for s in gpumd_input_atoms.get_chemical_symbols())


def test_load_gpumd_input():
    """Load all information from a GPUMD input file."""
    with open('xyz.in', 'w') as gpumd_input_f:
        gpumd_input_f.write(gpumd_input_text)

    species_types = {'Al': 0}
    gpumd_input_atoms, input_parameters, type_symbol_map =\
        load_xyz_input_gpumd('xyz.in', species_types=species_types)
    input_parameters_ref = {'N': 15, 'M': 1, 'cutoff': 0.1,
                            'use_triclinic': 0, 'has_velocity': 0,
                            'num_of_groups': 1}
    assert all(k in input_parameters for k in input_parameters_ref.keys())
    assert all(v == input_parameters[k] for k, v in
               input_parameters_ref.items())
    type_symbol_map_ref = {v: k for k, v in species_types.items()}
    assert all(k in type_symbol_map for k in type_symbol_map_ref.keys())
    assert all(v == type_symbol_map[k] for k, v in
               type_symbol_map_ref.items())


def test_gpumd_input_write():
    """Write a structure and read it back."""
    atoms = bulk('NiO', 'rocksalt', 4.813, cubic=True)

    # Test write and read
    atoms.write('xyz.in')
    readback = io.read('xyz.in')
    assert np.allclose(atoms.positions, readback.positions)
    assert np.allclose(atoms.cell, readback.cell)
    assert np.array_equal(atoms.numbers, readback.numbers)

    # Test write and load with groupings
    groupings = [[[i for i, s in
                   enumerate(atoms.get_chemical_symbols()) if s == 'Ni'],
                  [i for i, s in
                   enumerate(atoms.get_chemical_symbols()) if s == 'O']],
                 [[i] for i in range(len(atoms))]]
    groups = [[[j for j, group in enumerate(grouping) if i in group][0]
               for grouping in groupings] for i in range(len(atoms))]
    atoms.write('xyz.in', groupings=groupings)
    readback, input_parameters, _ = load_xyz_input_gpumd('xyz.in')
    assert input_parameters['num_of_groups'] == 2
    assert len(readback.info) == len(atoms)
    assert all(np.array_equal(
        readback.info[i]['groups'], np.array(groups[i])) for i in
        range(len(atoms)))

    # Test write and read with velocities
    velocities = np.array([[-0.3, 2.3, 0.7], [0.0, 0.3, 0.8],
                           [-0.6, 0.9, 0.1], [-1.7, -0.1, -0.5],
                           [-0.5, 0.0, 0.6], [-0.2, 0.1, 0.5],
                           [-0.1, 1.4, -1.9], [-1.0, -0.5, -1.2]])
    atoms.set_velocities(velocities)
    atoms.write('xyz.in')
    readback, input_parameters, _ = load_xyz_input_gpumd('xyz.in')
    assert input_parameters['has_velocity'] == 1
    assert np.allclose(readback.get_velocities(), atoms.get_velocities())
