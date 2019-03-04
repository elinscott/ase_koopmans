from __future__ import print_function
import os
from collections import OrderedDict
from collections import Mapping
from copy import deepcopy
from ase.atoms import Atoms
from ase.io.acemolecule import read_acemolecule_out
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.calculator import ReadError
from ase.calculators.calculator import Calculator
import numpy as np

class ACE(FileIOCalculator,Calculator):
    '''
    ACE-Molecule logfile reader
    '''
    name = 'ace'
    implemented_properties = ['energy', 'forces',
                              'geometry', 'excitation-energy']
#    system_changes = None
    # defaults is default value of ACE-input
    basic_list = [{
        'Type': 'Scaling', 'Scaling': '0.35', 'Basis': 'Sinc',
                  'Grid': 'Sphere',
                  'KineticMatrix': 'Finite_Difference', 'DerivativesOrder': '7',
                  'GeometryFilename': None, 'NumElectrons': None}
                  ]
    guess_list = []  # now not need this
    scf_list = [{
        'ExchangeCorrelation': {'XFunctional': 'GGA_X_PBE', 'CFunctional': 'GGA_C_PBE'},
        'NumberOfEigenvalues': None,
    }]

    force_list = [{'ForceDerivative': 'Potential'}]
    tddft_list = [{
        'SortOrbital': 'Order', 'MaximumOrder': '10',
        'ExchangeCorrelation': {'XFunctional': 'GGA_X_PBE', 'CFunctional': 'GGA_C_PBE'},\
    }]
    order_list = [0, 1, 2]
    cis_list = []
    cisd_list = []
    dda_list = []
    order_key_list = ['BasicInformation', 'Guess', 'Scf', 'Force', 'TDDFT','CIS', 'CISD', 'DDA', 'AVAS', 'LOCS', 'EnsembleKSScf', 'Charge', 'Analysis', 'DeltaScf', 'TAOScf', 'ExDFT', 'SigamScf']

    default_parameters = {'BasicInformation': basic_list, 'Guess': guess_list,
                          'Scf': scf_list, 'Force': force_list, 'TDDFT': tddft_list, 'order': order_list}
    parameters = default_parameters
    command = 'mpirun -np 1 ../ace PREFIX.inp > PREFIX.log'

    def __init__(
            self, restart=None, ignore_bad_restart_file=False,
            label='ace', atoms=None, command=None,
            basisfile=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command=command, **kwargs)

    def compare_parameters(self, parameters, val2):
        '''Replace parameters that users want'''
        for val_key, val_val in val2.items():
            for key, val in parameters.items():
                if val_key == key and (isinstance(val_val, str) or isinstance(val_val, float) or isinstance(val_val, int) or isinstance(val_val, list)):
                    parameters[key] = str(val2[key])
                elif (val_key == key and isinstance(val_val, dict)):
                    parameters[key] = self.compare_parameters(
                        parameters[key], val_val)

        return parameters

    def get_property(self, name, atoms=None, allow_calculation=True):
        '''Make input, xyz after that calculate and get_property(energy, forces, and so on)'''
        force_in_param = 0
        if(name=='forces'):
            if not 3 in self.parameters["order"]:
                self.parameters['order'].append(3)
                force_in_param = 1
                self.results = {}
                self.write_input(atoms)
        result = super().get_property(name, atoms, allow_calculation)
        if(force_in_param ==1):
            self.parameters['order'].pop()
           
        return result

    def set(self, **kwargs):
        changed_parameters = deepcopy(self.parameters)
        duplication = []
        if 'order' in kwargs:
            order_list = []
            for element in kwargs['order']:
                mod = 0
                if(element in self.order_key_list):
                    order_element = self.order_key_list.index(element)
                    order_list.append(order_element)
                    mod = 1
            if(mod == 1):
                kwargs['order'] = order_list
            changed_parameters['order'] = kwargs['order']
            for i in range(10):
                j = 0
                for value in kwargs['order']:
                    if(i in self.order_list and i == value):
                        j = j + 1
                        if(j > 1):
                            for num in range(j - 1):
                                changed_parameters[self.order_key_list[i]
                                                   ] += self.default_parameters[self.order_key_list[i]]

        for key in self.order_key_list:  # key : BasicInformation, Force, Scf and so on
            modified = False
            if key in kwargs.keys():  # kwargs.key() : BasicInformation, Force, ....
                if(isinstance(kwargs[key], dict)):
                    dict_to_list = []
                    dict_to_list.append(kwargs[key])
                    kwargs[key] = dict_to_list
                i = 0
                # kwargs[key] : basic_list, force_lsit ....
                for val in kwargs[key]:
                    element = self.compare_parameters(
                        changed_parameters[key][i], val)
                    if(element == changed_parameters[key][i]):
                        changed_parameters[key][i].update(val)
                    else:
                        duplication.append(element)
                        modified = True
                    i = i + 1
                    if(modified):
                        changed_parameters[key] = duplication
                        duplication = []
        self.parameters = changed_parameters
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
        self.parameters["BasicInformation"][0]["GeometryFilename"] = "{}_opt.xyz".format(
            self.label)
        self.parameters["BasicInformation"][0]["GeometryFormat"] = "xyz"
        self.write_acemolecule_input(inputfile, self.parameters)

        inputfile.close()

    def read_results(self):
        '''Read results from logfile '''
        filename = self.label + '.log'
        f = open(filename, "r")
        tddft = len(f.read().split("TDDFT"))
        if(tddft > 2):
            quantities = ['excitation-energy']
        else:
            quantities = ['energy', 'forces', 'atoms', 'excitation-energy']
        for value in quantities:
            self.results[value] = read_acemolecule_out(
                filename, quantity=value)

    def write_acemolecule_section(self, fpt, section, indent=0):
        for key, val in section.items():
            if(isinstance(val, str) or isinstance(val, int) or isinstance(val, float)):
                fpt.write('\t' * indent + str(key) + " " + str(val) + "\n")
            elif isinstance(val, dict):
                fpt.write('\t' * indent + "%% " + str(key) + "\n")
                indent = indent + 1
                self.write_acemolecule_section(fpt, val, indent)
                indent = indent - 1
                fpt.write('\t' * indent + "%% End\n")

        return

    def write_acemolecule_input(self, fpt, param2, indent=0):
        prefix = "    " * indent
        param = deepcopy(param2)
        for i in param['order']:
            fpt.write(prefix + "%% " + self.order_key_list[i] + "\n")
            section_list = param[self.order_key_list[i]]
            if(len(section_list) > 0):
                section = section_list.pop(0)
                self.write_acemolecule_section(fpt, section, 1)
            fpt.write("%% End\n")
        return


if __name__ == "__main__":
    ace = ACE()
    ace.set(ACEtemplate="test.a")
    ace.write_input(Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.0]]))
