.. module:: ase.calculators.ace_cal

========
ACE-Molecule
========

`ACE-Molecule <https://gitlab.com/aceteam.kaist/ACE-Molecule/wikis/home>`_  ACE-Molecule, whitch
stands for Advanced Computational Engine for Molecule, is a quantum chemistry package based on a 
real-space numerical grid. The package aims to carry out accurate and fast electronic structure 
calculations for large molecular systems. The present main features include ground-state DFT 
calculations using LDA, GGA, and hybrid functionals with an exact-exchange potential under the KLI 
approximation, atomic force calculations, and excited state calculations using 
linear-response TD-DFT, CIS, and CISD methods.

The ASE calculator is an interface to the ``ace`` executable.

Setup
=====
A simple calculation can be set up::

    import sys
    import os
    from ase.atoms import Atoms
    from ase.io import read
    from ase.calculators.ace_cal import ACE
    
    label = sys.argv[1]    
    mol= read('H2.xyz')
    ace = ACE(label=label)
    mol.set_calculator(ace)
    print (mol.get_potential_energy())

A Force calculation can be set up::
    
    import sys
    import os
    from ase.atoms import Atoms
    from ase.io import read
    from ase.calculators.ace_cal import ACE
    
    basic_list = [{'Pseudopotential' : {'Pseudopotential' : 1, 'Format' : 'upf', 'PSFilePath' : '/PATH/TO/UPF/FILES', 'PSFileSuffix' : '.pbe-theos.UPF'} } ]
    label = sys.argv[1]    
    mol= read('H2.xyz')
    order_list = [0, 1, 2, 3]
    ace = ACE(label=label, BasicInformation = basic_list, Order = order_list)
    mol.set_calculator(ace)
    print (mol.get_forces())
    

A Geometry optimization calculation can be set up:: 
    import sys
    import os
    from ase.atoms import Atoms
    from ase.io import read
    from ase.calculators.ace_cal import ACE
    from ase.optimize import BFGS

    basic_list = [{'Pseudopotential' : {'Pseudopotential' : 1, 'Format' : 'upf', 'PSFilePath' : '/PATH/TO/UPF/FILES', 'PSFileSuffix' : '.pbe-theos.UPF'} } ]
    label = sys.argv[1]    
    mol= read('H2.xyz')
    order_list = [0, 1, 2, 3]
    ace = ACE(label=label, BasicInformation = basic_list, Order = order_list)
    mol.set_calculator(ace)
    g_opt = BFGS(mol)
    g_opt.run(fmax=0.05)
    print ("OPT is end")

A TDDFT calculation can be set up ::
   import sys
   import os
   from ase.atoms import Atoms
   from ase.io import read
   from ase.calculators.ace_cal import ACE
   from ase.optimize import BFGS
   
   label = sys.argv[1]    
   mol= read('H2.xyz')
   order_list = [0, 1, 2, 4]
   scf_list = [dict(FinalDiagonalize = dict(NumberOfEigenvalues= 12))]
   ace = ACE(label=label, Scf= scf_list, Order = order_list)
   mol.set_calculator(ace)
   print (ace.get_property('excitation-energy', mol))
    
Parameters
==========

The calculator will interpret any of the documented options for ``ace``:
https://gitlab.com/aceteam.kaist/ACE-Molecule/tree/release-1.0.0/vars

default parameters exist. But If you want to obtain result of force 
calculation or TDDFT, you have to revise parameters. 
Parameters can be given as keywords and the calculator will put them into
the correct section of the input file.

The example of updating parameters::

    basic = [dict(Cell = 5.0, VerboseLevel = 2)]
    ace.set(BasicInformation = basic)


