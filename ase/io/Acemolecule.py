from __future__ import print_function
import numpy as np
import ase.units
from ase.atoms import Atoms
from ase.io.Acemoleculereader import Acemoleculereader as ACE_reader


def check_filename(filename):
    '''Adjust filename to match form'''
    check = str(filename).split("'")
    if(len(check) > 2):
        filename = check[1]
    return str(filename)

    
def read_Acemolecule_out(filename, quantity='atoms'):
    '''interface to Acemoleculereader and returns various quantities'''
    filename = check_filename(filename)
    data = ACE_reader(filename).data
    atom_symbol = np.array(data["Atomic_numbers"])
    positions = np.array(data["Positions"])
    atoms = Atoms(atom_symbol, positions=positions)

    f = open(filename, 'r')
    lines = f.read()
    energy_list = lines.split("Total energy ")
    energy_line = energy_list[len(energy_list) - 1]
    energy = float(energy_line.split('\n')[0].split('=')[1])
    geometry = zip(atom_symbol, positions)
        
    if(quantity == 'energy'):
        f.close()
        energy *= ase.units.Hartree
        # energy must be modified, hartree to eV
        return energy
    if(quantity == 'forces'):
        try:
            forces_lines = lines.split("total force in atomic unit.")[1]
            forces_line = forces_lines.split("======")[0].split('\n')
            forces = list()
            for i in range(2, len(forces_line) - 1):
                forces += [[float(forces_line[i].split()[3]),
                            float(forces_line[i].split()[4]),
                            float(forces_line[i].split()[5])]]
            convert = ase.units.Hartree / ase.units.Bohr
            forces = np.array(forces) * convert
        except:
            forces = None
        f.close()
        return forces
    if(quantity == 'geometry'):
        f.close()
        return geometry
    if(quantity == 'atoms'):
        f.close()
        return atoms


def read_Acemolecule_input(Label):
    '''Reads a Acemolecule input file'''
    filename = check_filename(Label)
    inputtemplate = open(filename, 'r')
    geometryfile_line = inputtemplate.read().split('GeometryFilename')[1]
    geometryfile = geometryfile_line.split('\n')[0].split()[0]
    xyzfile = open(geometryfile, 'r')
    atom_num = int(xyzfile.readline())
    atom_info = xyzfile.read().split('\n')[1:]
    atom_symbols = str()
    positions = []
    for i in range(atom_num):
        atom_symbols += atom_info[i].split()[0]
        x = atom_info[i].split()[1]
        y = atom_info[i].split()[2]
        z = atom_info[i].split()[3]
        positions.append((x, y, z))
    atoms = Atoms(atom_symbols, positions=positions)
    inputtemplate.close()
    xyzfile.close()
    return atoms
