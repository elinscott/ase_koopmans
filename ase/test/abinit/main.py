from ase.utils import workdir
from ase.build import bulk, molecule
from ase.calculators.abinit import Abinit


def run_test(atoms, name):
    dirname = 'test-abinit/{}'.format(name)
    with workdir(dirname, mkdir=True):
        header = 'test {} in {}'.format(name, dirname)
        print(header)
        print('=' * len(header))
        atoms.get_potential_energy()
        atoms.get_forces()
        eig = atoms.calc.results['eigenvalues']
        fermi = atoms.calc.results['fermilevel']
        print('fermi', fermi)
        print('eigshape', eig.shape)
    print()
    return atoms.calc.results


def abinit(**kwargs):
    kw = dict(ecut=150,
              chksymbreak=0,
              toldfe=1e-3)
    kw.update(kwargs)
    print(kw)
    return Abinit(**kw)


def test_si():
    atoms = bulk('Si')
    atoms.calc = abinit(nbands=4 * len(atoms))
    run_test(atoms, 'bulk')


def test_fe():
    atoms = bulk('Fe')
    atoms.set_initial_magnetic_moments([1])
    calc = abinit(nbands=8,
                  spinmagntarget=2.3,  # We'd rather have abinit find this number
                  kpts=[2, 2, 2])
    atoms.calc = calc
    run_test(atoms, 'bulk-spin')


def test_h2o():
    atoms = molecule('H2O', vacuum=2.5)
    atoms.calc = abinit(nbands=8)
    run_test(atoms, 'molecule')


def test_o2():
    atoms = molecule('O2', vacuum=2.5)
    atoms.calc = abinit(nbands=8,
                        spinmagntarget=2.0)
    run_test(atoms, 'molecule-spin')


test_si()
test_fe()
test_h2o()
test_o2()
