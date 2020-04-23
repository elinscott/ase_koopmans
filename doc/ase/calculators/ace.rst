.. module:: ase.calculators.acemolecule

============
ACE-Molecule
============

`ACE-Molecule <https://gitlab.com/aceteam.kaist/ACE-Molecule/wikis/home>`_ ACE-Molecule, whitch
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
    from ase.io import read
    from ase.calculators.acemolecule import ACE
    
    label = sys.argv[1]    
    mol= read('H2.xyz')
    basic_list = {'Cell' : 12.0}
    ace = ACE(label=label, BasicInformation = basic_list)
    mol.calc = ace
    print (mol.get_potential_energy())

A Force calculation can be set up::
    
    import sys
    from ase.io import read
    from ase.calculators.acemolecule import ACE
    
    basic_list = {'Cell' : 12.0 ,'Pseudopotential' : {'Pseudopotential' : 1, 'Format' : 'upf', 'PSFilePath' : '/PATH/TO/UPF/FILES', 'PSFileSuffix' : '.pbe-theos.UPF'} }
    label = sys.argv[1]    
    mol= read('H2.xyz')
    order_list = ["BasicInformation", "Guess", "Scf", "Force"]
    ace = ACE(label=label, BasicInformation = basic_list, order = order_list)
    mol.calc = ace
    print (mol.get_forces())
    

A Geometry optimization calculation can be set up:: 

    import sys
    from ase.io import read
    from ase.calculators.acemolecule import ACE
    from ase.optimize import BFGS

    basic_list = {'Cell' : 12.0, 'Pseudopotential' : {'Pseudopotential' : 1, 'Format' : 'upf', 'PSFilePath' : '/PATH/TO/UPF/FILES', 'PSFileSuffix' : '.pbe-theos.UPF'} }
    label = sys.argv[1]    
    mol= read('H2.xyz')
    order_list = ["BasicInformation", "Guess", "Scf", "Force"]
    ace = ACE(label=label, BasicInformation = basic_list, order = order_list)
    mol.calc = ace
    g_opt = BFGS(mol)
    g_opt.run(fmax=0.05)
    print ("OPT is end")

A TDDFT calculation can be set up::

   import sys
   from ase.io import read
   from ase.calculators.acemolecule import ACE
   from ase.optimize import BFGS

   label = sys.argv[1]
   mol= read('H2.xyz')
   basic_list = {'Cell' : 12.0}
   order_list = ["BasicInformation", "Guess", "Scf", "TDDFT"]
   scf_list = [dict(FinalDiagonalize = dict(NumberOfEigenvalues= 12))]
   ace = ACE(label=label, BasicInformation = basic_list, Scf= scf_list, order = order_list)
   mol.calc = ace
   print (ace.get_property('excitation-energy', mol))


Parameters
==========

The calculator will interpret any of the documented options for ``ace``:
https://gitlab.com/aceteam.kaist/ACE-Molecule/tree/master/vars

By default, this calculator sets simple SCF calculations (Including BasicInformation, Guess, and Scf section).
If you want to do force or excited state calculations, or change exchange-correlation functionals, you need to change input parameters.
Parameters can be given as keywords and the calculator will put them into the corresponding section of the input file.
Each (sub)section is represented as dictionary.

An example for updating parameters::

    basic = dict(Cell = 12.0, VerboseLevel = 2)
    ace.set(BasicInformation = basic)

An example for updating subsection parameters::

    ace.set(Scf = {"ExchangeCorrelation": {"XFunctional": "LDA_X", "CFunctional": "LDA_C_PZ"}})

