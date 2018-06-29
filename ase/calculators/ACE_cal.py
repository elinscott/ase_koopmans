from __future__ import print_function
import os
from collections import OrderedDict
from collections import Mapping
from copy import deepcopy
from re import sub
from ase.atoms import Atoms
from ase.io.Acemolecule import read_Acemolecule_out
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import ReadError


def update_nested(d, u):
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = update_nested(d.get(k, {}), v)
        else:
            d[k] = v
    return d


class OrderedParameters(OrderedDict):
    '''Dictionary for parameters.
    Since ACE-Molecule input section depends on order, we need OrderedDict.
    Nothing different with ase.calculators.calculator.Parameters
    '''
    def __init__(self, *args, **kwargs):
        super(OrderedParameters, self).__init__(*args, **kwargs)

    def __getattr__(self, key):
        if not key.startswith('_'):
            if key not in self:
                return OrderedDict.__getattribute__(self, key)
            return self[key]
        else:
            return super(OrderedParameters, self).__getattr__(key)

    def __setattr__(self, key, value):
        if not key.startswith('_'):
            self[key] = value
        else:
            return super(OrderedParameters, self).__setattr__(key, value)

    @classmethod
    def read(cls, filename):
        '''Read parameters from file.'''
        file = open(os.path.expanduser(filename))
        parameters = cls(eval(file.read()))
        file.close()
        return parameters

    def tostring(self):
        keys = sorted(self)
        return 'OrderedDict(' + ',\n     '.join(
            '{}={!r}'.format(key, self[key]) for key in keys) + ')\n'

    def write(self, filename):
        file = open(filename, 'w')
        file.write(self.tostring())
        file.close()


class ACE(FileIOCalculator):
    '''
    ACE-Molecule logfile reader
    '''
    name = 'ACE'
    implemented_properties = ['energy', 'forces', 'geometry']
    system_changes = None
    command = 'mpirun -np 1 ace PREFIX_opt.inp > PREFIX_opt.log'

    def __init__(
            self, restart=None, ignore_bad_restart_file=False,
            label='ACE', atoms=None, command=None,
            basisfile=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command=command, **kwargs)

    def get_default_parameters(self):
        pass

    def set(self, **kwargs):
        self.parameters = OrderedParameters()
        changed_parameters2 = OrderedParameters()
        if 'ACEtemplate' in kwargs:
            filename = kwargs.pop("ACEtemplate")
            changed_parameters2 = self.read_Acemolecule_inp(filename)

        changed_parameters = FileIOCalculator.set(self, **kwargs)
        update_nested(changed_parameters2, changed_parameters)
        self.parameters = changed_parameters2

        if changed_parameters:
            self.reset()
        return changed_parameters

    def read(self, label):
        FileIOCalculator.read(self, label)
        filename = self.label + "_opt.log"
        if not os.path.isfile(filename):
            raise ReadError
        self.read_results()

    def make_xyz_file(self, atoms):
        atoms.write(self.label + "_opt.xyz")

    def write_input(self, atoms, properties=None, system_changes=None):
        '''Writes the input file'''
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        inputfile = open(self.label + '_opt.inp', 'w')
        atoms.write("{}_opt.xyz".format(self.label))
        self.parameters["BasicInformation"]["GeometryFilename"] = "{}_opt.xyz".format(self.label)
        self.parameters["BasicInformation"]["GeometryFormat"] = "xyz"
        new_parameters = OrderedParameters()
        new_parameters["BasicInformation"] = OrderedDict()
        new_parameters["Guess"] = OrderedDict()
        update_nested(new_parameters, self.parameters)
        self.parameters = new_parameters
        self.write_Acemolecule_input(inputfile, self.parameters)

        inputfile.close()

    def read_results(self):
        filename = self.label + '_opt.log'
        quantities = ['energy', 'forces', 'atoms']
        for value in quantities:
            self.results[value] = read_Acemolecule_out(filename, quantity=value)

    def read_Acemolecule_inp(self, filename):
        with open(filename) as fpt:
            return self.parse_Acemolecule_inp(fpt)

    def parse_Acemolecule_inp(self, fpt, preserve_comment=False, sublist_name="root"):
        param = OrderedParameters()
        for line2 in fpt:
            line = line2.strip()
            line = line.replace('%%', '%% ')
            line = sub(r"(#+)", r"\1 ", line)
            if not preserve_comment:
                line = line.split('#')[0].strip()
            if len(line) == 0:
                continue
            L_line = line.split()

            if len(L_line) < 2:
                if not L_line[0] in param:
                    param[L_line[0]] = ""
            else:
                if L_line[0] == "%%":
                    if L_line[1] == "End":
                        return param
                    else:
                        param[L_line[1]] = self.parse_Acemolecule_inp(fpt, preserve_comment, L_line[1])
                elif L_line[0] in param:
                    if not type(param[L_line[0]]) is list:
                        param[L_line[0]] = [param[L_line[0]]]
                    param[L_line[0]].append(" ".join(L_line[1:]))
                else:
                    param[L_line[0]] = " ".join(L_line[1:])
        if not "root" == sublist_name:
            print(param)
            raise ReadError("Not matching ending block in " + sublist_name)
        return param

    def append_args(self, old_param, new_args):
        param = deepcopy(old_param)
        for key, val in new_args.items():
            key_list = key.split('.')
            param2 = param
        for key2 in key_list[:-1]:
            param2 = param2.setdefault(key2, OrderedDict())
        if val is not None:
            if type(val) is list:
                param2[key_list[-1]] = val
            else:
                param2[key_list[-1]] = str(val)
        return param

    def write_Acemolecule_input(self, fpt, param, indent=0):
        prefix = "    " * indent
        for key, val in param.items():
            if isinstance(val, Mapping):
                fpt.write(prefix + "%% " + str(key) + "\n")
                self.write_Acemolecule_input(fpt, param[key], indent + 1)
                fpt.write(prefix + "%% End\n")
            elif isinstance(val, list):
                for item in val:
                    if '#' in key:
                        fpt.write(prefix + str(key) + " " + str(item) + "\n")
                    else:
                        fpt.write(prefix + str(key) + "\t\t" + str(item) + "\n")
            else:
                if '#' in key:
                    fpt.write(prefix + str(key) + " " + str(val) + "\n")
                else:
                    fpt.write(prefix + str(key) + "\t\t" + str(val) + "\n")
        return

system_changes = None

if __name__ == "__main__":
    ace = ACE()
    ace.set(ACEtemplate="test.a")
    ace.write_input(Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.0]]))
