from ase_koopmans.io.opls import OPLSStructure
from ase_koopmans.visualize import view

s = OPLSStructure('172_mod.xyz')  # 172_mod.xyz if the file name for the structure above
view(s)  # view with real elements
elements = {'CT': 'Si', 'HC': 'H', 'H1': 'He'}
view(s.colored(elements))  # view with fake elements
