from gpaw import GPAW
calc = GPAW('bulk.Ag.gpw')
atoms = calc.get_atoms()
path = atoms.cell.bandpath('WLGXWK', density=10)
path.write('path.json')

calc.set(kpts=path, fixdensity=True, symmetry='off')

atoms.get_potential_energy()
bs = calc.band_structure()
bs.write('bs.json')
