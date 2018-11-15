from __future__ import division
import numpy as np

import ase.units as units
from ase.calculators.calculator import Calculator, all_changes
from ase.data import atomic_masses

# Electrostatic constant
k_c = units.Hartree * units.Bohr

# Force field parameters
q_me = 0.206
q_c = 0.247
q_n = -0.453
sigma_me = 3.775
sigma_c = 3.650
sigma_n = 3.200
sigma = np.array([sigma_me, sigma_c, sigma_n])
epsilon_me = 0.7824 * units.kJ / units.mol
epsilon_c = 0.544 * units.kJ / units.mol
epsilon_n = 0.6276 * units.kJ / units.mol
epsilon = np.array([epsilon_me, epsilon_c, epsilon_n])
r_mec = 1.458
r_cn = 1.157
r_men = r_mec + r_cn

# Variables needed to distribute forces from central C to Methyl and N
c_me = r_cn / r_men
c_n = r_mec / r_men
m_me = atomic_masses[6] + 3 * atomic_masses[1]
m_c = atomic_masses[6]
m_n = atomic_masses[7]
m_cn = m_c * m_n
m_mec = m_c * m_me
m_men = m_n * m_me
n_me = c_me / (c_me**2 * m_cn + c_n**2 * m_mec + m_men)
n_n = c_n / (c_me**2 * m_cn + c_n**2 * m_mec + m_men)


def wrap(D, cell, pbc):
    """Wrap distances to nearest neighbor (minimum image convention)."""
    shift = np.zeros_like(D)
    for i, periodic in enumerate(pbc):
        if periodic:
            d = D[:, i]
            L = cell[i]
            shift[:, i] = (d + L / 2) % L - L / 2 - d
    return shift


def combine_lj_lorenz_berthelot(sigma, epsilon):
    """Combine LJ parameters according to the Lorenz-Berthelot rule"""
    sigma_c = np.zeros((len(sigma), len(sigma)))
    epsilon_c = np.zeros_like(sigma_c)

    for ii in range(len(sigma)):
        sigma_c[:, ii] = (sigma[ii] + sigma) / 2
        epsilon_c[:, ii] = (epsilon[ii] * epsilon) ** 0.5
    return sigma_c, epsilon_c


class ACN(Calculator):
    implemented_properties = ['energy', 'forces']
    nolabel = True
    pcpot = None

    def __init__(self, rc=5.0, width=1.0):
        """Three-site potential for acetonitrile.
           
        Parameters from: 
        http://dx.doi.org/10.1080/08927020108024509

        Correct atom sequence should be:
        MeCNMeCN ... MeCN or NCMeNCMe ... NCMe 

        Forces are redistributed from the central C atom to N 
        and Methyl according to a scheme by Ciccotti et al. 
        (http://dx.doi.org/10.1080/00268978200100942) for MD
        propagation of rigid linear triatomic molecules. Use 
        class FixBondLengthsLinear in constraints.py to do MD 
        with this calculator.

        rc: float
            Cutoff radius for Coulomb part.
        width: float
            Width for cutoff function for Coulomb part.
        """
        self.rc = rc
        self.width = width
        self.forces = None
        Calculator.__init__(self)

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        R = self.atoms.positions.reshape((-1, 3, 3))
        pbc = self.atoms.pbc
        cell = self.atoms.cell.diagonal()
        nm = len(R)

        assert (self.atoms.cell == np.diag(cell)).all(), 'not orthorhombic'
        assert ((cell >= 2 * self.rc) | ~pbc).all(), 'cutoff too large'

        charges = self.get_virtual_charges(atoms[:3])

        energy = 0.0
        self.forces = np.zeros((3 * nm, 3))

        for m in range(nm - 1):
            Dmm = R[m + 1:, 1] - R[m, 1]
            # MIC PBCs
            shift = wrap(Dmm, cell, pbc)
        
            # Smooth cutoff
            Dmm += shift
            d2 = (Dmm**2).sum(1)
            d = d2**0.5
            cut, dcut = self.cutoff(d)

            # LJ parameters
            sigma_c, epsilon_c = combine_lj_lorenz_berthelot(sigma, epsilon)

            for j in range(3):
                D = R[m + 1:] - R[m, j] + shift[:, np.newaxis]
                r2 = (D**2).sum(axis=2)
                r = r2**0.5
                # Coulomb interactions
                e = charges[j] * charges / r * k_c
                energy += np.dot(cut, e).sum()
                F = (e / r2 * cut[:, np.newaxis])[:, :, np.newaxis] * D
                Fmm = -(e.sum(1) * dcut / d)[:, np.newaxis] * Dmm
                self.forces[(m + 1) * 3:] += F.reshape((-1, 3))
                self.forces[m * 3 + j] -= F.sum(axis=0).sum(axis=0)
                self.forces[(m + 1) * 3 + 1::3] += Fmm
                self.forces[m * 3 + 1] -= Fmm.sum(0)
                # LJ interactions
                c6 = (sigma_c[:, j]**2 / r2)**3
                c12 = c6**2
                e = 4 * epsilon_c[:, j] * (c12 - c6)
                energy += np.dot(cut, e).sum()
                F = (24 * epsilon_c[:, j] * (2 * c12 -
                     c6) / r2 * cut[:, np.newaxis])[:, :, np.newaxis] * D
                Fmm = -(e.sum(1) * dcut / d)[:, np.newaxis] * Dmm
                self.forces[(m + 1) * 3:] += F.reshape((-1, 3))
                self.forces[m * 3 + j] -= F.sum(axis=0).sum(axis=0)
                self.forces[(m + 1) * 3 + 1::3] += Fmm
                self.forces[m * 3 + 1] -= Fmm.sum(0)
                
        if self.pcpot:
            e, f = self.pcpot.calculate(np.tile(charges, nm),
                                        self.atoms.positions)
            energy += e
            self.forces += f
        
        fr = self.redistribute_forces(self.forces)
        self.forces = fr

        self.results['energy'] = energy
        self.results['forces'] = self.forces
 
    def redistribute_forces(self, fo):
        fr = np.zeros_like(fo)
        Z = self.atoms.numbers
        if Z[0] == 7:
            n = 0
            me = 2
        else:
            n = 2
            me = 0
        assert (Z[n::3] == 7).all(), 'Incorrect atoms sequence'
        assert (Z[1::3] == 6).all(), 'Incorrect atoms sequence'
        
        # N
        fr[n::3, :] = ((1 - n_n * m_mec * c_n) * fo[n::3, :] - 
                       n_n * m_cn * c_me * fo[me::3, :] + 
                       n_n * m_men * fo[1::3, :])
        # Me 
        fr[me::3, :] = ((1 - n_me * m_cn * c_me) * fo[me::3, :] - 
                        n_me * m_mec * c_n * fo[n::3, :] + 
                        n_me * m_men * fo[1::3, :])
 
        return fr                 

    def get_molcoms(self, nm):      
        molcoms = np.zeros((nm, 3))      
        for m in range(nm): 
            molcoms[m] = self.atoms[m * 3:(m + 1) * 3].get_center_of_mass() 
        return molcoms
 
    def cutoff(self, d): 
        x1 = d > self.rc - self.width  
        x2 = d < self.rc 
        x12 = np.logical_and(x1, x2)    
        y = (d[x12] - self.rc + self.width) / self.width  
        cut = np.zeros(len(d))  # cutoff function    
        cut[x2] = 1.0     
        cut[x12] -= y**2 * (3.0 - 2.0 * y)    
        dtdd = np.zeros(len(d))     
        dtdd[x12] -= 6.0 / self.width * y * (1.0 - y)        
        return cut, dtdd

    def embed(self, charges):
        """Embed atoms in point-charges."""
        self.pcpot = PointChargePotential(charges)
        return self.pcpot

    def check_state(self, atoms, tol=1e-15):
        system_changes = Calculator.check_state(self, atoms, tol)
        if self.pcpot and self.pcpot.mmpositions is not None:
            system_changes.append('positions')
        return system_changes

    def add_virtual_sites(self, positions):
        return positions  # no virtual sites

    def get_virtual_charges(self, atoms):
        charges = np.empty(len(atoms))
        Z = atoms.numbers
        if Z[0] == 7:
            n = 0
            me = 2
        else:
            n = 2
            me = 0
        assert (Z[n::3] == 7).all(), 'Incorrect atoms sequence'
        assert (Z[1::3] == 6).all(), 'Incorrect atoms sequence'
        charges[me::3] = q_me
        charges[1::3] = q_c
        charges[n::3] = q_n
        return charges
    

class PointChargePotential:
    def __init__(self, mmcharges):
        """Point-charge potential for ACN.

        Only used for testing QMMM.
        """
        self.mmcharges = mmcharges
        self.mmpositions = None
        self.mmforces = None

    def set_positions(self, mmpositions):
        self.mmpositions = mmpositions

    def calculate(self, qmcharges, qmpositions):
        energy = 0.0
        self.mmforces = np.zeros_like(self.mmpositions)
        qmforces = np.zeros_like(qmpositions)
        for C, R, F in zip(self.mmcharges, self.mmpositions, self.mmforces):
            d = qmpositions - R
            r2 = (d**2).sum(1)
            e = units.Hartree * units.Bohr * C * r2**-0.5 * qmcharges
            energy += e.sum()
            f = (e / r2)[:, np.newaxis] * d
            qmforces += f
            F -= f.sum(0)
        self.mmpositions = None
        return energy, qmforces

    def get_forces(self, calc):
        return self.mmforces
