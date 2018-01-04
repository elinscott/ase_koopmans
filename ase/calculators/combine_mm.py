from __future__ import print_function
import numpy as np
from ase.calculators.calculator import Calculator
from ase.calculators.qmmm import wrap
from ase import units
import copy

k_c = units.Hartree * units.Bohr
#k_c = 332.1 * units.kcal / units.mol

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
    
    Maybe it can combine n MM calculators in the future? """

    implemented_properties = ['energy', 'forces']
    def __init__(self, idx, apm1, apm2, calc1, calc2, vdw, rc=7.0, width=1.0):
        self.idx = idx
        self.apm1 = apm1  # atoms per mol
        self.apm2 = apm2
        ## self.molidx = molidx

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

        self.cell = atoms.cell
        self.pbc = atoms.pbc


    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.atoms1 is None:
            self.initialize(atoms)

        pos1 = atoms.positions[self.mask]
        pos2 = atoms.positions[~self.mask]
        self.atoms1.set_positions(pos1)
        self.atoms2.set_positions(pos2)

        # positions and charges for the coupling term, which should
        # include virtual charges and sites: 
        xpos1 = self.atoms1.calc.add_virtual_sites(pos1)
        xpos2 = self.atoms2.calc.add_virtual_sites(pos2)

        xc1 = self.atoms1.calc.get_virtual_charges(self.atoms1)
        xc2 = self.atoms2.calc.get_virtual_charges(self.atoms2)
        
        xpos1 = xpos1.reshape((-1, self.apm1, 3))
        xpos2 = xpos2.reshape((-1, self.apm2, 3)) 
        
        # shift for qmmm not used yet.  
        shift = np.array([0, 0, 0])

        e_c, f_c = self.coulomb(xpos1, xpos2, xc1, xc2, shift)

        # PBCs wrt total box should now also be applied to subsys 1
        # which is different from the qmmm method, for which the LJ was made.
        # so prewrap atoms1 here. 
        cell = atoms.cell.diagonal()
        pos = self.atoms1.get_positions()
        for i, periodic in enumerate(atoms.pbc):
            if periodic:
                d = pos[:, i]
                L = cell[i]
                d  = (d + L ) % L - L  

        watoms1 = self.atoms1.copy()
        watoms1.set_positions(pos)

        e_vdw, f1, f2 = self.vdw.calculate(watoms1, self.atoms2, shift)
        f_vdw = np.zeros((len(atoms), 3))
        f_vdw[self.mask] += f1
        f_vdw[~self.mask] += f2

        # internal energy, forces of each subsystem:
        f12 = np.zeros((len(atoms), 3))
        e1 = self.atoms1.get_potential_energy()
        fi1 = self.atoms1.get_forces()

        e2 = self.atoms2.get_potential_energy()
        fi2 = self.atoms2.get_forces()

        f12[self.mask] += fi1
        f12[~self.mask] += fi2

        self.results['energy'] = e_c + e_vdw + e1 + e2
        self.results['forces'] = f_c + f_vdw + f12

    def get_virtual_charges(self, atoms):
        vc = np.zeros(len(self.atoms))
        # this can break, IF there is virtual sites. XXX 
        vc1 = self.atoms1.calc.get_virtual_charges(atoms[self.mask])
        vc2 = self.atoms2.calc.get_virtual_charges(atoms[~self.mask])
        vc[self.mask] = vc1
        vc[~self.mask] = vc2

        return vc

    def add_virtual_sites(self, positions): 
        vs = np.zeros(len(self.atoms), 3)
        # this can break, IF there is virtual sites. XXX 
        vs1 = self.atoms1.calc.add_virtual_sites()
        vs2 = self.atoms2.calc.add_virtual_sites()
        vs[self.mask] = vs1
        vs[~self.mask] = vs2

        return vc

    def coulomb(self, xpos1, xpos2, xc1, xc2, shift):
        energy = 0.0
        forces = np.zeros((len(xc1)+len(xc2), 3))

        self.xpos1 = xpos1
        self.xpos2 = xpos2

        R1 = xpos1  
        R2 = xpos2
        F1 = np.zeros_like(R1)
        F2 = np.zeros_like(R2)
        C1 = xc1.reshape((-1, self.apm1))
        C2 = xc2.reshape((-1, self.apm2))
        # Vectorized evaluation is not as trivial when apm1 != apm2.
        # This is pretty inefficient, but for ~1-5 counter ions as region 1
        # it should not matter much ..
        # There is definetely room for improvements here.
        cell = self.cell.diagonal()
        for m1, (r1, c1) in enumerate(zip(R1, C1)):
            for m2, (r2, c2) in enumerate(zip(R2, C2)):
                r00 = r2[0] - r1[0]
                shift = np.zeros(3)
                for i, periodic in enumerate(self.pbc):
                    if periodic:
                        L = cell[i]
                        shift[i] = (r00[i] + L / 2.) % L - L / 2. - r00[i]
                r00 += shift  

                d00 = (r00**2).sum()**0.5
                t = 1
                dtdd = 0
                if d00 > self.rc:
                    continue 
                elif d00 > self.rc - self.width:
                    y = (d00 - self.rc + self.width) / self.width
                    t -= y**2 * (3.0 - 2.0 *y)  
                    dtdd = r00 * 6 * y * (1.0 - y) / (self.width * d00) 

                for a1 in range(self.apm1):
                    for a2 in range(self.apm2):
                        r = r2[a2] - r1[a1] + shift
                        d2 = (r**2).sum()
                        d = d2**0.5
                        e = k_c * c1[a1] * c2[a2] / d
                        energy += t * e

                        F1[m1, a1] -= t * (e / d2) * r 
                        F2[m2, a2] += t * (e / d2) * r

                        F1[m1, 0] -= dtdd * e  
                        F2[m2, 0]  += dtdd * e 


        F1 = F1.reshape((-1, 3))
        F2 = F2.reshape((-1, 3))

        # Redist forces but dont save forces in org calculators
        atoms1 = self.atoms1.copy()
        atoms1.calc = copy.copy(self.calc1)
        atoms1.calc.atoms = atoms1
        F1 = atoms1.calc.redistribute_forces(F1)
        atoms2 = self.atoms2.copy()
        atoms2.calc = copy.copy(self.calc2)
        atoms2.calc.atoms = atoms2
        F2 = atoms2.calc.redistribute_forces(F2)

        forces = np.zeros((len(self.atoms), 3))
        forces[self.mask] = F1
        forces[~self.mask] = F2

        return energy, forces



