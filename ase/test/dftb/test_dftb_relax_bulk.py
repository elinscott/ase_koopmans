def test_dftb_relax_bulk():
    import os
    from ase.test import require
    from ase.test.testsuite import datafiles_directory
    from ase.build import bulk
    from ase.calculators.dftb import Dftb
    from ase.optimize import QuasiNewton
    from ase.constraints import ExpCellFilter

    require('dftb')

    os.environ['DFTB_PREFIX'] = datafiles_directory

    calc = Dftb(label='dftb',
                kpts=(3,3,3),
                Hamiltonian_SCC='Yes')

    atoms = bulk('Si')
    atoms.set_calculator(calc)

    ecf = ExpCellFilter(atoms)
    dyn = QuasiNewton(ecf)
    dyn.run(fmax=0.01)

    e = atoms.get_potential_energy()
    assert abs(e - -73.150819) < 1., e
