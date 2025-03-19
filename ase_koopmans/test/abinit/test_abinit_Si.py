from ase_koopmans.build import bulk
from ase_koopmans.units import Ry


def test_abinit_Si(abinit_factory):
    atoms = bulk('Si')
    atoms.calc = abinit_factory.calc(
        label='Si',
        nbands=8,
        ecut=10 * Ry,
        kpts=[4, 4, 4],
        toldfe=1.0e-2,
    )
    e = atoms.get_potential_energy()
    print(e)
