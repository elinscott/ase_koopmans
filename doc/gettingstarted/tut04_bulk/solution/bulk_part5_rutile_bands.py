from gpaw import GPAW

calc = GPAW('groundstate.rutile.gpw')
atoms = calc.get_atoms()
path = atoms.cell.bandpath(density=7)
path.write('path.rutile.json')

calc.set(kpts=path, fixdensity=True,
         symmetry='off')

atoms.get_potential_energy()
bs = calc.band_structure()
bs.write('bs.rutile.json')
