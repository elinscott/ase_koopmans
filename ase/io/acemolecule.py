from __future__ import print_function
import numpy as np
import ase.units
from ase.atoms import Atoms
from ase.io.acemoleculereader import Acemoleculereader as ACE_reader



    
def read_acemolecule_out(filename, quantity='atoms'):
    '''interface to Acemoleculereader and returns various quantities'''
    data = ACE_reader(filename).data
    atom_symbol = np.array(data["Atomic_numbers"])
    positions = np.array(data["Positions"])
    atoms = Atoms(atom_symbol, positions=positions)

    f = open(filename, 'r')
    lines = f.readlines()
    geometry = zip(atom_symbol, positions)
        
    if(quantity == 'excitation-energy'):
        f.close()
        # ee is excitation-energy
        ee = 1
        return ee
    if(quantity == 'energy'):
        energy = 0
        for i in range(len(lines)-1,1,-1):
            line = lines[i].split()
            if(len(line)>2):
                if(line[0] == 'Total' and line[1] == 'energy'):
                    energy = float(line[3])
                    break
        f.close()
        energy *= ase.units.Hartree
        # energy must be modified, hartree to eV
        return energy
    if(quantity == 'forces'):
        for i in range(len(lines)-1,1,-1):
            if (lines[i] == '!================================================\n'):
                endline_num = i
            if (lines[i] == '! Atom           x         y         z\n'):
                forces = []
                startline_num = i
                for j in range(startline_num+1, endline_num):
                    forces += [[float(lines[j].split()[3]), 
                                float(lines[j].split()[4]), 
                                float(lines[j].split()[5])]]
                convert = ase.units.Hartree / ase.units.Bohr
                forces = np.array(forces) * convert
                break
            if(i == 30):
                forces = None
                break
        f.close()
        return forces
    if(quantity == 'geometry'):
        f.close()
        return geometry
    if(quantity == 'atoms'):
        f.close()
        return atoms


def read_acemolecule_input(Label):
    '''Reads a Acemolecule input file'''
    filename = Label
    inputtemplate = open(filename, 'r')
    lines = inputttemplate.readlines()
    for line in lines:
        if(len(line.split('GeometryFilename')) > 1):
            line.split()[2] = geometryfile
            break
    atoms = read(geometryfile, format='xyz')
    inputtemplate.close()
    return atoms
