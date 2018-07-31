from __future__ import print_function
import sys
from ase.data import chemical_symbols


def check_filename(filename):
    check = str(filename).split("'")
    if(len(check) > 2):
        filename = check[1]
    return str(filename)


class Acemoleculereader:
    
    def auto_type(self, data):
        ''' tries to determine type'''
        try:
            return float(data)
        except ValueError:
            pass

        try:
            ds = data.split(",")
            array = []

            for d in ds:
                array.append(float(d))

            return array
        except ValueError:
            pass

        return data

    def __init__(self, filename):
        filename = check_filename(filename)
        self.parse(filename)

    def parse(self, filename):
        ''' Read atoms geometry '''
        f = open(filename, 'r')
        lines = f.read()
        geometry_lines = lines.split("  Atoms  =====================")[-1]
        geo_lines = geometry_lines.split("===============")[0].split('\n')
        self.data = []
        atoms = []
        positions = []
        new_dict = {}
        for i in range(1, len(geo_lines) - 1):
            atomic_number = geo_lines[i].split()[0]
            atoms.append(str(chemical_symbols[int(atomic_number)]))
            x = geo_lines[i].split()[1]
            y = geo_lines[i].split()[2]
            z = geo_lines[i].split()[3]
            positions.append((x, y, z))
        new_dict["Atomic_numbers"] = atoms
        new_dict["Positions"] = positions
        self.data = new_dict


if __name__ == "__main__":
    from ase.io import read as ACE_read
    Label = str(sys.argv[1].split('.inp')[0])
    system_changes = None
    a = ACE_read(Label + '.inp', format='acemolecule-input')
    print(a.get_positions())
    
    filename = Label + '.log'
