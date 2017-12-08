from __future__ import print_function
import numpy as np
from ase.calculators.calculator import Calculator
from ase.calculators.qmmm import LJInteractionsGeneral as LJG
from ase import units

k_c = units.Hartree * units.Bohr

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
        ## self.molidx1 = self.molidx[self.mask]
        ## self.molidx2 = self.molidx[~self.mask]

        atoms.constraints = constraints

        self.atoms1.calc = self.calc1
        self.atoms2.calc = self.calc2


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
        
        # Do pbc stuff to get shift
        shift = np.array([0, 0, 0])

        e_c, f_c = self.coulomb(xpos1, xpos2, xc1, xc2, shift)

        e_vdw, f1, f2 = self.vdw.calculate(self.atoms1, self.atoms2, shift)
        f_vdw = np.zeros((len(atoms), 3))
        f_vdw[self.mask] = f1
        f_vdw[~self.mask] = f2

        # internal energy, forces of each subsystem:
        f12 = np.zeros((len(atoms), 3))
        e1 = self.atoms1.get_potential_energy()
        f1 = self.atoms1.get_forces()
        e2 = self.atoms2.get_potential_energy()
        f2 = self.atoms2.get_forces()
        f12[self.mask] = f1
        f12[~self.mask] = f2

        self.results['energy'] = e_c + e_vdw + e1 + e2
        self.results['forces'] = f_c + f_vdw + f12

    def get_virtual_charges(self, atoms):
        vc = np.zeros(len(self.atoms))
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
        '''  Not done... 
            '''
        e = 0.0
        forces = np.zeros((len(xc1)+len(xc2), 3))

        ## idx1 = np.array(self.molidx1)
        ## idx2 = np.array(self.molidx2)
        
        ## mols1_pos = [pos1[np.where(idx1==i)[0]] for i in set(idx1)]
        ## mols2_pos = [pos2[np.where(idx2==i)[0]] for i in set(idx2)]

        R1 = xpos1.reshape((-1, self.apm1, 3))  # molwise for cutoff
        R2 = xpos2.reshape((-1, self.apm2, 3))
        C1 = xc1.reshape((-1, self.apm1))
        C2 = xc2.reshape((-1, self.apm2))

        for m1, c1 in zip(R1, C1):  
            for m2, c2 in zip(R2, C2):
                d00 = (sum((m1[0] - m2[0])**2))**0.5
                if d00 > self.rc:  # molwise cutoff from 1st atom in each mol
                    continue
                t = 1
                if d00 > self.rc - self.width:
                    y = (d00 - self.rc + self.width) / self.width
                    t -= y**2 * (3.0 - 2.0 *y)  # same value for entire mols
                for a1 in range(self.apm1):
                    for a2 in range(self.apm2):
                        r = m1[a1] - m2[a2]
                        d2 = (r**2).sum()
                        d = d2**0.5 
                        e += t * k_c * c1[a1] * c2[a2] / d

#        for m in mols1_pos:
#            DOO = R[m + 1:, o] - R[m, o]
#            #shift = np.zeros_like(DOO)
#            #for i, periodic in enumerate(pbc):
#            #    if periodic:
#            #        L = cell[i]
#            #        shift[:, i] = (DOO[:, i] + L / 2) % L - L / 2 - DOO[:, i]
#            DOO += shift
#            d2 = (DOO**2).sum(1)
#            d = d2**0.5
#            x1 = d > self.rc - self.width
#            x2 = d < self.rc
#            x12 = np.logical_and(x1, x2)
#            y = (d[x12] - self.rc + self.width) / self.width
#            t = np.zeros(len(d))  # cutoff function
#            t[x2] = 1.0
#            t[x12] -= y**2 * (3.0 - 2.0 * y)
#            dtdd = np.zeros(len(d))
#            dtdd[x12] -= 6.0 / self.width * y * (1.0 - y)
#
#            for j in range(atoms_in_molecule):
#                D = R[m + 1:] - R[m, j] + shift[:, np.newaxis]
#                r2 = (D**2).sum(axis=2)
#                r = r2**0.5
#                e = charges[j] * charges / r * units.Hartree * units.Bohr
#                XXX
#                energy += np.dot(t, e).sum()
#                F = (e / r2 * t[:, np.newaxis])[:, :, np.newaxis] * D
#                FOO = -(e.sum(1) * dtdd / d)[:, np.newaxis] * DOO
#                forces[(m + 1) * 3 + o::3] += FOO
#                forces[m * 3 + o] -= FOO.sum(0)
#                forces[(m + 1) * 3:] += F.reshape((-1, 3))
#                forces[m * 3 + j] -= F.sum(axis=0).sum(axis=0)
#       
#

        # ima need some force redistribution here... 
        return e, np.zeros((len(self.atoms),3))


    def coulomb_nocutoff(self, c, shift):
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

        return energy, ft

    def coulomb_nomolidx(self, c, shift):
        '''  Not done... 
            '''
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
            r = r2**0.5

            # cutoff
            x1 = r > self.rc - self.width
            x2 = r < self.rc
            x12 = np.logical_and(x1, x2)
            y = (r[x12] - self.rc + self.width) / self.width
            t = np.zeros(len(r))  # cutoff function
            t[x2] = 1.0
            t[x12] -= y**2 * (3.0 - 2.0 * y)
            dtdr = np.zeros(len(r))
            dtdr[x12] -= 6.0 / self.width * y * (1.0 - y)

            XXX

            e = units.Hartree * units.Bohr * C * c1 * r**-1 
            energy += np.dot(t, e).sum()
            f = (e / r2)[:, np.newaxis] * d
            f1 += f
            F -= f.sum(0)
       
        ft[self.mask] = f1
        ft[~self.mask] = f2

        return energy, ft


