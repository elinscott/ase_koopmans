def test_orca():
    from ase_koopmans.optimize import BFGS
    from ase_koopmans.atoms import Atoms
    from ase_koopmans.calculators.orca import ORCA


    atoms = Atoms('OHH',
                  positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])

    atoms.calc = ORCA(label='water',
                      orcasimpleinput='BLYP def2-SVP')


    opt = BFGS(atoms)
    opt.run(fmax=0.05)

    final_energy = atoms.get_potential_energy()
    print(final_energy)

    assert abs(final_energy + 2077.24420) < 1.0
