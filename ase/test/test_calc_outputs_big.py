import pytest
import numpy as np

from ase.outputs import Properties, all_properties


@pytest.fixture
def rng():
    return np.random.RandomState(17)


def test_properties_big(rng):
    nspins, nkpts, nbands = 2, 3, 5
    natoms = 4

    results = dict(
        natoms=natoms,
        energy=rng.rand(),
        free_energy=rng.rand(),
        energies=rng.rand(natoms),
        forces=rng.rand(natoms, 3),
        stress=rng.rand(6),
        stresses=rng.rand(natoms, 6),
        nspins=nspins,
        nkpts=nkpts,
        nbands=nbands,
        eigenvalues=rng.rand(nspins, nkpts, nbands),
        occupations=rng.rand(nspins, nkpts, nbands),
        fermi_level=rng.rand(),
        ibz_kpoints=rng.rand(nkpts, 3),
        kpoint_weights=rng.rand(nkpts),
    )

    props = Properties()
    for name, value in results.items():
        props.setvalue(name, value)
    #assert set(out) == all_properties, all_properties ^ set(out)

    for name in all_properties:
        assert name in props, name
        obj = props[name]
        #obj = getattr(out, name)
        print(name, obj)
