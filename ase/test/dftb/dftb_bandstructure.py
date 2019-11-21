from ase.calculators.dftb import Dftb
from ase.build import bulk

calc = Dftb(label='dftb',
            kpts=(3,3,3),
            Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance=1e-5,
            Hamiltonian_MaxAngularMomentum_Si='d')

atoms = bulk('Si')
atoms.set_calculator(calc)
atoms.get_potential_energy()

efermi = calc.get_fermi_level()
assert abs(efermi - -2.90086680996455) < 1.

# DOS does not currently work because of 
# missing "get_k_point_weights" function
#from ase.dft.dos import DOS
#dos = DOS(calc, width=0.2)
#d = dos.get_dos()
#e = dos.get_energies()
#print(d, e)

calc = Dftb(atoms=atoms,
            label='dftb',
            kpts={'path':'WGXWLG', 'npoints':50},
            Hamiltonian_SCC='Yes',
            Hamiltonian_MaxSCCIterations=1,
            Hamiltonian_ReadInitialCharges='Yes',
            Hamiltonian_MaxAngularMomentum_Si='d')

atoms.set_calculator(calc)
calc.calculate(atoms)

calc.results['fermi_levels'] = [efermi]
bs = calc.band_structure()
bs.plot(filename='bandstructure.png')
