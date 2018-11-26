from __future__ import print_function
import os
from collections import OrderedDict
from collections import Mapping
from copy import deepcopy
from re import sub
from ase.atoms import Atoms
from ase.io.acemolecule import read_acemolecule_out
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import ReadError
import numpy as np


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
    implemented_properties = ['energy', 'forces', 'geometry','excitation-energy']
    system_changes = None
    parameters = OrderedParameters()
    user_parameters = OrderedParameters()
    # defaults is default value of ACE-input
    defaults = OrderedParameters([('BasicInformation',\
                                 OrderedParameters([('Label', '/home/khs/work/asetest/H2test'),\
                                                     ('Type', 'Scaling'), ('Scaling', '1.0'), ('Basis', 'Sinc'),\
                                                     ('Cell', '3.0'), ('Grid', 'Sphere'),\
                                                     ('KineticMatrix', 'Finite_Difference'), ('DerivativesOrder', '7'),\
                                                     ('GeometryFormat', 'xyz'), ('SpinMultiplicity', '1.0'), ('Centered', '0'),\
                                                     ('OccupationMethod', 'ZeroTemp'), ('Polarize', '0'),\
                                                     ('Pseudopotential',\
                                                     OrderedParameters([('Pseudopotential', '1'),\
                                                                       ('Format', 'upf'),\
                                                                       ('PSFilePath', '/home/khs/DATA/UPF')\
                                                                       ,('PSFileSuffix', '.pbe-theos.UPF')\
                                                                       ])),\
                                                     ('GeometryFilename', 'test_opt.xyz'), ('NumElectrons', '2')])),\
                                  ('Guess',\
                                 OrderedParameters([('InitialGuess', '1'), ('InitialFilePath', '/home/khs/DATA/UPF'),\
                                                    ('InitialFileSuffix', '.pbe-theos.UPF')])),\
                                  ('Scf',\
                                 OrderedParameters([('IterateMaxCycle', '150'), ('ConvergenceType', 'Energy'),\
                                                    ('ConvergenceTolerance', '0.00001'), ('EnergyDecomposition', '1'),\
                                                    ('ComputeInitialEnergy', '1'),\
                                                    ('Diagonalize',\
                                                   OrderedParameters([('DiagonalizeMaxIter', '10'), ('Tolerance', '0.000001'),\
                                                                      ('FullOrthogonalize', '1')])),\
                                                    ('ExchangeCorrelation',\
                                                   OrderedParameters([('CalGradientUsingDensity', '1'),\
                                                                      ('XFunctional', 'GGA_X_PBE'), ('CFunctional', 'GGA_C_PBE')\
                                                                      ])),\
                                                    ('Mixing',\
                                                   OrderedParameters([('MixingMethod', '1'), ('MixingType', 'Density'),\
                                                                      ('MixingParameter', '0.5'), ('PulayMixingParameter', '0.1')\
                                                                      ])),\
                                                    ('SolvationModel',\
                                                   OrderedParameters([('Solvent', 'water'), ('Area', '1.0'),\
                                                    ('SolverType', 'None')])),\
                                                    ('NumberOfEigenvalues', '4')])),\
                                  ('Force',\
                                 OrderedParameters([('ForceDerivative', 'Potential')]))])
    command = 'mpirun -np 1 ace PREFIX.inp > PREFIX.log'

    def __init__(
            self, restart=None, ignore_bad_restart_file=False,
            label='ACE', atoms=None, command=None,
            basisfile=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command=command, **kwargs)
        
    def compare_parameters(self,parameters, key2, val2):    
        '''Replace and Add parameters that users want'''
        compare = 0
#        print("???????????????????????????")
#        print (parameters)
        for key, val in parameters.items():
            if(key == key2):
                parameters[key] = val2
                print("a parameter is updated")
                print (type([1,1]))
                return 10
            if( compare != 10 and type(parameters[key]) != type('1') and type(parameters[key]) != type([1,1]) \
                and type(parameters[key]) != type(None) ):
                ##### type is different depanding on python version.
#                print (type(parameters[key]))
                compare = self.compare_parameters(parameters[key],key2, val2)
        return compare
    def get_property(self, name, atoms=None,  allow_calculation=True):
        '''Make input, xyz after that calculate and get_property(energy, forces, and so on)'''
        
        if name not in self.implemented_properties:
            print ('{} property not implemented'.format(name))
        if atoms is None:
            atoms = self.atoms
            system_changes = []
        else:
            system_changes = self.check_state(atoms)
            if system_changes:
                self.reset()
        if (name == 'excitation-energy'):
            self.parameters =\
            OrderedParameters([('BasicInformation',\
                              OrderedParameters([('Type', 'Scaling'), ('Scaling', '1.0'), ('Basis', 'Sinc'), ('Grid', 'Atoms'),\
                                                 ('Radius', ['4.0', '4.0']), ('AbsoluteRadius', '1'),\
                                                 ('GeometryFilename', 'xyz/benzene.xyz'), ('GeometryFormat', 'xyz'),\
                                                 ('SpinMultiplicity', '1.0'), ('Centered', '1'), ('OccupationMethod', 'ZeroTemp'),\
                                                 ('Polarize', '0'), ('NumElectrons', '30'), ('KineticMatrix', 'Finite_Difference'),\
                                                 ('DerivativesOrder', '9'), ('ShallowCopyOrbitals', '1'),('ParallelGridAtoms', '1'),\
                                                 ('Pseudopotential',\
                                                OrderedParameters([('Pseudopotential', '1'), ('Format', 'upf'),\
                                                                   ('UsingDoubleGrid', '0'), ('NonlocalRmax', '2.0'),\
                                                                   ('Nonlocalthreshold', '0.000001'), \
                                                                   ('PSFilePath', '/home/khs/DATA/UPF/'), \
                                                                   ('PSFileSuffix', '.pbe-theos.UPF')]))])),\
                               ('Guess',\
                              OrderedParameters([('InitialGuess', '3'), ('InitialFilePath', 'pbe.gs/'), ('Info', 'benzene_opt.log'),\
                                                 ('InfoFiletype', 'ACE'), ('NumberOfEigenvalues', '40')])),\
                               ('TDDFT',\
                              OrderedParameters([('SortOrbital', 'Order'), ('MaximumOrder', '30'), ('Gradient', '1'),\
                                                 ('ExchangeCorrelation', \
                                                OrderedParameters([('CalGradientUsingDensity', '0'), ('XFunctional', '101'),\
                                                                   ('CFunctional', '130')])),\
                                                 ('OrbitalInfo',\
                                                OrderedParameters([('ExchangeCorrelation', \
                                                                  OrderedParameters([('CalGradientUsingDensity', '0'),\
                                                                                     ('XFunctional', '101'),\
                                                                                     ('CFunctional', '130')]))])),\
                                                 ('xckernelwithoutEXX',\
                                                OrderedParameters([('ExchangeCorrelation',\
                                                                 OrderedParameters([('CalGradientUsingDensity', '0')]))]))]))])
            
        if 'modified' in self.parameters:
            del self.parameters['modified'] 
        for key, val in self.user_parameters.items(): 
            none_value = self.compare_parameters(self.parameters, key, val)   
            if(none_value == 0):
                self.parameters[key] = val
        if 'system_changes' in self.parameters:
            del self.parameters['system_changes']
        self.write_input(atoms) 


        if name not in self.results:
            if not allow_calculation:
                return None
            self.calculate(atoms, [name], system_changes)
        
        if name not in self.results:
            # For some reason the calculator was not able to do what we want,
            # and that is OK.
            print ('{} not present in this calculation'.format(name))
        result = self.results[name]
        if isinstance(result, np.ndarray):
            result = result.copy()
        
        return result

    def set(self, **kwargs):
        self.parameters = OrderedParameters()
        changed_parameters2 = self.defaults
        if 'ACEtemplate' in kwargs:
            filename = kwargs.pop("ACEtemplate")
            changed_parameters2 = self.read_acemolecule_inp(filename)
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        update_nested(changed_parameters2, changed_parameters)
        self.parameters = changed_parameters2
        if 'modified' in self.parameters:
            self.user_parameters = kwargs.pop('modified')
            del self.parameters['modified'] 
        if changed_parameters:
            self.reset()
        return changed_parameters

    def read(self, label):
        FileIOCalculator.read(self, label)
        filename = self.label + ".log"
        if not os.path.isfile(filename):
            raise ReadError
        self.read_results()

    def make_xyz_file(self, atoms):
        atoms.write("{}_opt.xyz".format(self.label))

    def write_input(self, atoms, properties=None, system_changes=None):
        '''Writes the input file and xyz file'''
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        inputfile = open(self.label + '.inp', 'w')
        self.make_xyz_file(atoms)
        self.parameters["BasicInformation"]["GeometryFilename"] = "{}_opt.xyz".format(self.label)
        self.parameters["BasicInformation"]["GeometryFormat"] = "xyz"
        new_parameters = OrderedParameters()
        new_parameters["BasicInformation"] = OrderedDict()
        new_parameters["Guess"] = OrderedDict()
        update_nested(new_parameters, self.parameters)
        self.parameters = new_parameters
        self.write_acemolecule_input(inputfile, self.parameters)

        inputfile.close()

    def read_results(self):
        '''Read results from logfile '''
        filename = self.label + '.log'
        f= open(filename,"r")
        tddft = len(f.read().split("TDDFT"))
        if(tddft>2):
            quantities = ['excitation-energy']
        else:
            quantities = ['energy','forces', 'atoms','excitation-energy']
        for value in quantities:
            self.results[value] = read_acemolecule_out(filename, quantity=value)

    def read_acemolecule_inp(self, filename):
        with open(filename) as fpt:
            return self.parse_acemolecule_inp(fpt)

    def parse_acemolecule_inp(self, fpt, preserve_comment=False, sublist_name="root"):
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
                        param[L_line[1]] = self.parse_acemolecule_inp(fpt, preserve_comment, L_line[1])
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

    def write_acemolecule_input(self, fpt, param, indent=0):
        prefix = "    " * indent
        for key, val in param.items():
            if isinstance(val, Mapping):
                fpt.write(prefix + "%% " + str(key) + "\n")
                self.write_acemolecule_input(fpt, param[key], indent + 1)
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
