import numpy as np

from ase.outputs import CalculatorOutputs, all_properties

rng = np.random.RandomState(17)

nspins, nkpts, nbands = 2, 3, 5
natoms = 4

results = dict(
    natoms=natoms,
    energy=rng.rand(),
    free_energy=rng.rand(),
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

for name in all_properties:
    obj = getattr(out, name)
    print(name, obj)
