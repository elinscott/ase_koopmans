def test_Al_dmol():
    from ase_koopmans.build import bulk
    from ase_koopmans.calculators.dmol import DMol3

    atoms = bulk('Al')
    calc = DMol3()
    atoms.calc = calc
    atoms.get_potential_energy()
    atoms.get_forces()
