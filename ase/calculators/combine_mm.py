from __future__ import print_function
import numpy as np
from ase.calculators.calculator import Calculator
from ase.calculators.qmmm import LJInteractionsGeneral as LJG
from ase import units

class CombineMM(Calculator):
    """Hopefully a calculator that combines two MM calculators 
    (TIPnP, ACN, counterions). 

    It needs to: 

    - Define what parts belong to what calc
    - Do PBC stuff and remember cutoffs
    - Obtain forces and energies from both subsets
    - Calculate forces and energies from their interaction:
        - electrostatic
        - vdw
    - Return values
    - Be embeddable with the Embedding class, so it needs:
        - get_virtual_charges()
        - some way around virtual_molecule_size
        - make molecule_size into a list, and ask if its
            longer than 1, then do smarter reshapes to get
            the correct positions arrays. 

    
    Maybe it can combine n MM calculators in the future? """

    implemented_properties = ['energy', 'forces']
    def __init__(self, idx, calc1, calc2, vdw, rc=7.0, width=1.0):
        self.idx = idx
        self.rc = rc 
        self.width = width

        self.atoms1 = None
        self.atoms2 = None
        self.mask = None

        self.calc1 = calc1
        self.calc2 = calc2

        self.vdw = vdw

        Calculator.__init__(self)

    def initialize(self, atoms):
        self.mask = np.zeros(len(atoms), bool)
        self.mask[self.idx] = True

        constraints = atoms.constraints
        atoms.constraints = []
        self.atoms1 = atoms[self.mask]
        self.atoms2 = atoms[~self.mask]
        atoms.constraints = constraints

        self.atoms1.calc = self.calc1
        self.atoms2.calc = self.calc2

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)


        if self.atoms1 is None:
            self.initialize(atoms)

        self.atoms1.set_positions(atoms.positions[self.mask])
        self.atoms2.set_positions(atoms.positions[~self.mask])
        
        #### Do pbc stuff to get shift
        shift = np.array([0, 0, 0])

        c = self.get_virtual_charges(atoms)
        e_c, f_c = self.coulomb(c, shift) 

        e_vdw, f1, f2 = self.vdw.calculate(self.atoms1, self.atoms2, shift)
        f_vdw = np.zeros((len(atoms), 3))
        f_vdw[self.mask] = f1
        f_vdw[~self.mask] = f2

        self.results['energy'] = e_c + e_vdw
        self.results['forces'] = f_vdw

    def get_virtual_charges(self, atoms):
        vc = np.zeros(len(self.atoms))
        vc1 = self.atoms1.calc.get_virtual_charges(atoms[self.mask])
        vc2 = self.atoms2.calc.get_virtual_charges(atoms[~self.mask])
        vc[self.mask] = vc1
        vc[~self.mask] = vc2

        return vc

    def add_virtual_sites(self, positions): 
        vs = np.zeros(len(self.atoms), 3)
        vs1 = self.atoms1.calc.add_virtual_sites()
        vs2 = self.atoms2.calc.add_virtual_sites()
        vs[self.mask] = vs1
        vs[~self.mask] = vs2

        return vc

    def coulomb(self, c, shift):
        ''' Needs cutoff and new variable names. 
            Maybe look into tipnp instead'''
        energy = 0.0
        f1 = np.zeros((len(self.atoms1), 3))
        f2 = np.zeros((len(self.atoms2), 3))
        ft = np.zeros((len(self.atoms), 3))

        c1 = c[self.mask]
        c2 = c[~self.mask]

        pos1 = self.atoms1.positions
        pos2 = self.atoms2.positions

        for C, R, F in zip(c2, pos2, f2):
            d = pos1 - R
            r2 = (d**2).sum(1)
            e = units.Hartree * units.Bohr * C * r2**-0.5 * c1
            energy += e.sum()
            f = (e / r2)[:, np.newaxis] * d
            f1 += f
            F -= f.sum(0)
#        self.mmpositions = None  # this one is weird
       
        ft[self.mask] = f1
        ft[~self.mask] = f2

        return energy, f

#    def vdw(self, atoms1, atoms2, shift):
#        lj = LJG(sig1, eps1, sig2, eps2)
#        e, f1, f2 = lj.calculate(atoms1, atoms2, shift)
#        f = np.zeros((len(self.atoms), 3))
#        f[self.mask] = f1
#        f[~self.mask] = f2
#        return e, f2



