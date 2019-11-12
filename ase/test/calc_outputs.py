import numpy as np

from ase.outputs import CalculatorOutputs, all_properties

rng = np.random.RandomState(17)

nspins, nkpts, nbands = 2, 3, 5
natoms = 4

results = dict(
    natoms=natoms,
    energy=rng.rand(),
    free_energy=rng.rand(),
    energies=rng.rand(natoms),
    forces=rng.rand(natoms, 3),
    stress=rng.rand(6),
    nspins=nspins,
    nkpts=nkpts,
    nbands=nbands,
    eigenvalues=rng.rand(nspins, nkpts, nbands),
    occupations=rng.rand(nspins, nkpts, nbands),
    fermi_level=rng.rand(),
    ibz_kpoints=rng.rand(nkpts, 3),
    kpoint_weights=rng.rand(nkpts),
)

out = CalculatorOutputs(results)
#assert set(out) == all_properties, all_properties ^ set(out)

for name in all_properties:
    assert name in out, name
    obj = getattr(out, name)
    print(name, obj)
