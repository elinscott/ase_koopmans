from __future__ import print_function
import os
from collections import OrderedDict
from collections import Mapping
from copy import deepcopy
#from re import sub
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
    name = 'ace'
    implemented_properties = ['energy', 'forces', 'geometry','excitation-energy']
    system_changes = None
    user_parameters = OrderedParameters()
    # defaults is default value of ACE-input
    basic_list =[{
                  'Type' : 'Scaling', 'Scaling' : '0.35', 'Basis' : 'Sinc',\
                  'Cell' : '3.0', 'Grid' : 'Sphere',\
                  'KineticMatrix': 'Finite_Difference', 'DerivativesOrder' : '7',\
                  'GeometryFilename': None, 'NumElectrons': None}
                 ]
    guess_list = [] #now not need this
    scf_list = [ {
                 'ExchangeCorrelation' : {'XFunctional' : 'GGA_X_PBE', 'CFunctional' : 'GGA_C_PBE'} ,\
                 'NumberOfEigenvalues': None,
                 }]
    
    force_list = [{ 'ForceDerivative' : 'Potential'   }  ]
#    solvation_list = [{'SolvationLibrary' : 'PCMSolver', 'Solvent' : 'water', 'Area' : '0.2', 'SolverType' : 'CPCM' }] 
    tddft_list = [{
    'SortOrbital': 'Order', 'MaximumOrder' : '10',\
            'ExchangeCorrelation' :  {'XFunctional' : 'GGA_X_PBE', 'CFunctional' : 'GGA_C_PBE'},\
#            'OrbitalInfo' : {'ExchangeCorrelation' : {'XFunctional' : 'GGA_X_PBE', 'CFunctional' : 'GGA_C_PBE'}},\
            }] 
    order_list = [0,1,2]
    order_key_list = ['BasicInformation', 'Guess', 'Scf','Force', 'TDDFT' ]


    default_parameters = {'BasicInformation': basic_list, 'Guess' : guess_list, 'Scf':scf_list, 'Force' : force_list, 'TDDFT': tddft_list ,'Order' : order_list} 
    parameters = default_parameters
    command = 'mpirun -np 1 ../ace PREFIX.inp > PREFIX.log'

    def __init__(
            self, restart=None, ignore_bad_restart_file=False,
            label='ace', atoms=None, command=None,
            basisfile=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command=command, **kwargs)
        
    def compare_parameters(self,parameters, key2, val2):    
        '''Replace parameters that users want'''
        for val_key, val_val in val2.items():
            for key, val in parameters.items():
                if val_key==key and (isinstance(val_val,str) or isinstance(val_val, float) or isinstance(val_val, int) or isinstance(val_val, list) ):
                    parameters[key] = str(val2[key])
                elif (val_key==key and isinstance(val_val,dict)):
                    parameters[key] = self.compare_parameters(parameters[key], key, val_val)            
        
        
        return parameters


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
            
#        if 'modified' in self.parameters:
#            del self.parameters['modified'] 
#        for key, val in self.user_parameters.items(): 
#            none_value = self.compare_parameters(self.parameters, key, val)   
#            if(none_value == 0):
#                self.parameters[key] = val
#        if 'system_changes' in self.parameters:
#            del self.parameters['system_changes']
        print("write input run")
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
#        self.parameters = OrderedParameters()
        changed_parameters = deepcopy(self.parameters)
#        if 'ACEtemplate' in kwargs:
#            filename = kwargs.pop("ACEtemplate")
#            changed_parameters2 = self.read_acemolecule_inp(filename)
#        changed_parameters = FileIOCalculator.set(self, **kwargs)
#        update_nested(changed_parameters2, changed_parameters)
#        self.parameters = changed_parameters2
        duplication = []
        if 'Order' in kwargs:
            changed_parameters['Order'] = kwargs['Order']
            for i in range(10):
                j = 0
                for value in kwargs['Order']:
                    if(i in self.order_list and i==value):
                        j= j+1
                        if(j>1):
                            for num in range(j-1):
                                changed_parameters[self.order_key_list[i]] += self.default_parameters[self.order_key_list[i]]
                            
        for key in self.order_key_list: #### key : BasicInformation, Force, Scf and so on
            modified = False
            print(kwargs)
            if key  in kwargs.keys(): ##### kwargs.key() : In Basic, Cell, GeometryFilename, ....
                i=0
                for val in kwargs[key]: ########## kwargs[key] : basic_list, force_lsit ....
                    element = self.compare_parameters(changed_parameters[key][i], key, val)
                    if(element == self.parameters[key][i]):
#                        print("yes")
#                        print(element)
                        print(self.parameters[key][i])
#                        print("yes end")
                        changed_parameters[key][i].update(val) 
                    else:
#                        print("duplication")
                        duplication.append(element)
                        modified = True
                    i= i+1
            if(modified):
                changed_parameters[key] = duplication
#        print("in_set")
#        print(changed_parameters)
        self.parameters = changed_parameters
#        if changed_parameters:
#            self.reset()
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
        print("write_input")
        print(self.parameters["BasicInformation"])
        print("BasciInfor end")
        self.parameters["BasicInformation"][0]["GeometryFilename"] = "{}_opt.xyz".format(self.label)
        self.parameters["BasicInformation"][0]["GeometryFormat"] = "xyz"        
        #new_parameters = OrderedParameters()
        #new_parameters["BasicInformation"] = OrderedDict()
        #new_parameters["Guess"] = OrderedDict()
        #update_nested(new_parameters, self.parameters)
        #self.parameters = new_parameters
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

#    def read_acemolecule_inp(self, filename):
#        with open(filename) as fpt:
#            return self.parse_acemolecule_inp(fpt)

#    def parse_acemolecule_inp(self, fpt, preserve_comment=False, sublist_name="root"):
#        param = OrderedParameters()
#        for line2 in fpt:
#            line = line2.strip()
#            line = line.replace('%%', '%% ')
#            line = sub(r"(#+)", r"\1 ", line)
#            if not preserve_comment:
#                line = line.split('#')[0].strip()
#            if len(line) == 0:
#                continue
#            L_line = line.split()
#
#            if len(L_line) < 2:
#                if not L_line[0] in param:
#                    param[L_line[0]] = ""
#            else:
#                if L_line[0] == "%%":
#                    if L_line[1] == "End":
#                        return param
#                    else:
#                        param[L_line[1]] = self.parse_acemolecule_inp(fpt, preserve_comment, L_line[1])
#                elif L_line[0] in param:
#                    if not type(param[L_line[0]]) is list:
#                        param[L_line[0]] = [param[L_line[0]]]
#                    param[L_line[0]].append(" ".join(L_line[1:]))
#                else:
#                    param[L_line[0]] = " ".join(L_line[1:])
#        if not "root" == sublist_name:
#            raise ReadError("Not matching ending block in " + sublist_name)
#        return param

#    def append_args(self, old_param, new_args):
#        param = deepcopy(old_param)
#        for key, val in new_args.items():
#            key_list = key.split('.')
#            param2 = param
#        for key2 in key_list[:-1]:
#            param2 = param2.setdefault(key2, OrderedDict())
#        if val is not None:
#            if type(val) is list:
#                param2[key_list[-1]] = val
#            else:
#                param2[key_list[-1]] = str(val)
#        return param

    def write_acemolecule_section(self, fpt, section,indent = 0):
        for key, val in section.items():
            if(isinstance(val,str) or isinstance(val,int) or isinstance(val,float)):
                fpt.write('\t'*indent + str(key) + " " + str(val) + "\n")
            elif isinstance(val, dict):
                fpt.write('\t'*indent +"%% " + str(key) + "\n")
                indent = indent+1
                self.write_acemolecule_section(fpt,val,indent)
                indent = indent-1
                fpt.write('\t'*indent +"%% End\n")
   
        return

    def write_acemolecule_input(self, fpt, param2, indent=0):
        prefix = "    " * indent
        param = deepcopy(param2)
        for i in param['Order']:
            fpt.write(prefix+"%%"+self.order_key_list[i] +"\n")
            section_list = param[self.order_key_list[i]]
            if(len(section_list) > 0):
                section = section_list.pop(0)
                print("section start")
                print(section)
                print("section end")
                self.write_acemolecule_section(fpt,section,1)
            fpt.write("%% End\n") 
        return


if __name__ == "__main__":
    ace = ACE()
    ace.set(ACEtemplate="test.a")
    ace.write_input(Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.0]]))
