def test_ace():
    import os
    import unittest

    from ase import Atoms
    from ase.calculators.acemolecule import ACE

    label = "test"
    mol = Atoms('H2',[(0, 0, 0),(0, 0, 0.7)])
    basic = [dict(Cell= '5.0')]
    if "ASE_ACE_COMMAND" not in os.environ:
        raise unittest.SkipTest('$ASE_ACE_COMMAND not defined') 
    ace = ACE(label=label, BasicInformation=basic)
    mol.calc = ace
    mol.get_forces()
    #mol.get_potential_energy()
