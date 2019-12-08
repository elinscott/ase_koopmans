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
    run_test(atoms, 'bulk-si')


def test_au():
    atoms = bulk('Au')
    atoms.calc = abinit(nbands=10 * len(atoms))
    run_test(atoms, 'bulk-au')


def _test_fe_anyocc(name, **kwargs):
    atoms = bulk('Fe')
    atoms.set_initial_magnetic_moments([1])
    calc = abinit(nbands=8,
                  kpts=[2, 2, 2], **kwargs)
    atoms.calc = calc
    run_test(atoms, name)

def test_fe_fixed_magmom():
    _test_fe_anyocc('bulk-spin-fixmagmom', spinmagntarget=2.3)

def test_fe_any_magmom():
    _test_fe_anyocc('bulk-spin-anymagmom', occopt=7)

def test_h2o():
    atoms = molecule('H2O', vacuum=2.5)
    atoms.calc = abinit(nbands=8)
    run_test(atoms, 'molecule')


def test_o2():
    atoms = molecule('O2', vacuum=2.5)
    atoms.calc = abinit(nbands=8,
                        spinmagntarget=2.0)
    run_test(atoms, 'molecule-spin')


def test_big():
    atoms = bulk('Au') * (2, 2, 2)
    atoms.rattle(stdev=0.01)
    atoms.symbols[:2] = 'Cu'
    atoms.calc = abinit(nbands=len(atoms) * 7,
                        kpts=[8, 8, 8])
    run_test(atoms, 'big')

def test_many():
    atoms = bulk('Ne', cubic=True) * (4, 2, 2)
    atoms.rattle(stdev=0.01)
    atoms.calc = abinit(nbands=len(atoms) * 5)
    run_test(atoms, 'manyatoms')

#test_many()
#test_big()
test_si()
test_au()
test_fe_fixed_magmom()
test_fe_any_magmom()
test_h2o()
test_o2()
