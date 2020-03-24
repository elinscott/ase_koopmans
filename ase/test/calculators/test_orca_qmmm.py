def test_orca_qmmm():
    from ase.calculators.tip4p import TIP4P, epsilon0, sigma0
    from ase.calculators.orca import ORCA
    from ase.calculators.qmmm import EIQMMM, LJInteractions
    from ase.data import s22

    atoms = s22.create_s22_system('Water_dimer')


    qmcalc = ORCA(label='water',
                  orcasimpleinput='BLYP def2-SVP')

    lj = LJInteractions({('O', 'O'): (epsilon0, sigma0)})


    atoms.calc = EIQMMM(selection=[0, 1, 2],
                        qmcalc=qmcalc,
                        mmcalc=TIP4P(),
                        interaction=lj,
                        output='orca_qmmm.log')


    e = atoms.get_potential_energy()

    assert abs(e + 2077.45445852) < 1.0
