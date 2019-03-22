import numpy as np
import ase.units
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from ase.data import chemical_symbols


def parse_geometry(self, filename):
    '''Read atoms geometry from ACE-Molecule log file and put it to self.data.
    Parameters
    ==========
    filename: ACE-Molecule log file.

    Returns
    =======
    Dictionary of parsed geometry data.
    retval["Atomic_numbers"]: list of atomic numbers
    retval["Positions"]: list of [x, y, z] coordinates for each atoms.
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
        start_line = 0
        end_line = 0
        for i in range(len(lines)):
            if lines[i] == '====================  Atoms  =====================\n':
                start_line = i
            if start_line != '' and len(lines[i].split('=')) == '3':
                end_line = i
                break
        atoms = []
        positions = []
        for i in range(start_line + 1, end_line):
            atomic_number = lines[i].split()[0]
            atoms.append(str(chemical_symbols[int(atomic_number)]))
            xyz = [float(n) for n in lines[i].split()[1:4]]
            positions.append(xyz)

        return {"Atomic_numbers": atoms, "Positions": positions}


def read_acemolecule_out(filename, quantity='atoms'):
    '''Interface to ACEMoleculeReader and return values for corresponding quantity
    Parameters
    ==========
    filename: ACE-Molecule log file.
    quantity: One of atoms, energy, forces, excitation-energy.

    Returns
    =======
     - quantity = 'excitation-energy':
       returns None. This is placeholder function to run TDDFT calculations without IndexError.
     - quantity = 'energy':
       returns energy as float value.
     - quantity = 'forces':
       returns force of each atoms as numpy array of shape (natoms, 3).
     - quantity = 'atoms':
       returns ASE atoms object.
    '''
    data = parse_geometry(filename)
    atom_symbol = np.array(data["Atomic_numbers"])
    positions = np.array(data["Positions"])
    atoms = Atoms(atom_symbol, positions=positions)

    with open(filename, 'r'):
        lines = f.readlines()
    energy = 0
    forces = None
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    atoms.set_calculator(calc)

    if quantity == 'excitation-energy':
        return None

    if quantity == 'energy':
        for i in range(len(lines) - 1, 1, -1):
            line = lines[i].split()
            if len(line) > 2:
                if line[0] == 'Total' and line[1] == 'energy':
                    energy = float(line[3])
                    break
        energy *= ase.units.Hartree
        # energy must be modified, hartree to eV
        return energy

    if quantity == 'forces':
        for i in range(len(lines) - 1, 1, -1):
            if ('!============================' in lines[i]):
                endline_num = i
            if ('! Atom        ' in lines[i]):
                forces = []
                startline_num = i
                for j in range(startline_num + 1, endline_num):
                    forces += [[float(lines[j].split()[3]),
                                float(lines[j].split()[4]),
                                float(lines[j].split()[5])]]
                convert = ase.units.Hartree / ase.units.Bohr
                forces = np.array(forces) * convert
                break
            if i == 30:
                forces = None
                break
        return forces

    if quantity == 'geometry':
        geometry = zip(atom_symbol, positions)
        return geometry

    if quantity == 'atoms':
        return atoms


def read_acemolecule_input(label):
    '''Reads a ACE-Molecule input file'''
    filename = label
    inputtemplate = open(filename, 'r')
    lines = inputtemplate.readlines()
    inputtemplate.close()
    for line in lines:
        if len(line.split('GeometryFilename')) > 1:
            geometryfile = line.split()[2]
            break
    atoms = read(geometryfile, format='xyz')
    return atoms


if __name__ == "__main__":
    import sys
    from ase.io import read as ACE_read
    Label = str(sys.argv[1].split('.inp')[0])
    system_changes = None
    a = ACE_read(Label + '.inp', format='acemolecule-input')

    filename = Label + '.log'
