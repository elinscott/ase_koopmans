from __future__ import print_function
import numpy as np
import sys
import os
from collections import Mapping
from copy import deepcopy
from re import sub
from ase.atoms import Atoms
from ase.io import read
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.ace_cal import ACE
from ase.calculators.ace_cal import OrderedParameters
from collections import OrderedDict
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
Label = str(sys.argv[1].split('.inp')[0])
system_changes=None
print ("Fine)")
filename=Label+'.inp'
################################# Command is "python ACE_opt.py [file_name(absoulte_path)]" or "qsub ase_jobscript()"



############################### user have to change this part ####################################
# SP : ground state Single Point E calculation, OPT : geomery optimization, TDDFT : TDDFT calculation when existing Cubefile
# SP+TDDFT : If you don't have Cubefile, makes Cube file. Afther than, TDDFT calculation start.  
mode = 'SP' 

mol= read('/home/khs/work/asetest/xyz/H2.xyz') # xyz must be writed as absoulte_path and
                                               # also all other values are writed as absoulte_paht
order_list = [0,1,2]
ace = ACE(label=Label,command = 'mpirun -np 4 /home/khs/hs_file/programs/ACE-Molecule/ace PREFIX.inp > PREFIX.log') ## You have to revise path of executable file.

basic = [dict(Cell= '4.0',VerboseLevel=1)]
ace.set(BasicInformation=basic)
###############################About below part, user doen't change anything###########################################################################

mol.set_calculator(ace)

if (mode == "SP" or mode == 'OPT'):
    print ('Single Point Energy : ',  ace.get_property('energy', mol))
    print ("SP calculation is end")
    if(mode == 'OPT'):
        g_opt = BFGS(mol)
        g_opt.run(fmax=0.05)
        print ("OPT is end")


if (mode == 'TDDFT'):
    print (ace.get_property('excitation-energy', mol))
    print ("TDDFT calculation is end")


if(mode == 'SP+TDDFT'):  #### Maybe this part is doesn't work. I'm revising this part.
    print (ace.get_property('energy', mol))
    print ("Cube files are made")
    mol= read('/home/khs/hs_file/programs/ACE-Molecule/aseinput/tddft_example/xyz/benzene.xyz')
    order_list = [0,1,4]
    ace = ACE(label=Label,command = 'mpirun -np 4 /home/khs/hs_file/programs/ACE-Molecule/ace PREFIX.inp > PREFIX.log', Order = order_list)
    mol.set_calculator(ace)
    print (ace.get_property('excitation-energy', mol))
    print ("TDDFT calculation is end")





###################################################### these are not important, just scientific digit change to digit.
def e_to_value(value):
   state=str(value).split("e-")
   integer=list()
   check =str(value).split("e+")
   if (len(check)==2):
       return check[0]
   if (len(state[0].split('-'))==2):
       if(len(state)==2):
           if(int(state[1])>7):
               return '0.0'
           return_value='-0.'
           state[0]=state[0].split('-')[1]
           for i in range(1,int(state[1])):
               return_value +='0'
           integer = state[0].split('.')
#            for i in range (len(integer)):
#                return_value+=integer[i]
#            return return_value
       else:
           return value
   else:
        if(len(state)==2):
            if(int(state[1])>7):
                return '0.0'
            return_value='0.'
            for i in range(1,int(state[1])):
                return_value +='0'
            integer = state[0].split('.')
            for i in range( len(integer)):
                return_value+=integer[i]
            return return_value
        else:
            return value

if(mode == 'OPT'):
    print ('================================')
        #print e
    position_list= mol.get_positions()
    distance=[]
    for position in position_list:
        position =str(position).split('[')[1].split(']')[0].split()
        distance.append(float(e_to_value(position[0])))
        distance.append(float(e_to_value(position[1])))
        distance.append(float(e_to_value(position[2])))
        print ( position[0]+' '+position[1]+' '+position[2] )
    value = (distance[0]-distance[3])**2+(distance[1]-distance[4])**2+(distance[2]-distance[5])**2
    value2= value**0.5
    print ( value2 )
    print ('================================')

