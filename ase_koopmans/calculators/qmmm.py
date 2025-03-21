import numpy as np

from ase_koopmans.calculators.calculator import Calculator
from ase_koopmans.data import atomic_numbers
from ase_koopmans.utils import convert_string_to_fd


class SimpleQMMM(Calculator):
    """Simple QMMM calculator."""

    implemented_properties = ['energy', 'forces']

    def __init__(self, selection, qmcalc, mmcalc1, mmcalc2, vacuum=None):
        """SimpleQMMM object.

        The energy is calculated as::

                    _          _          _
            E = E  (R  ) - E  (R  ) + E  (R   )
                 QM  QM     MM  QM     MM  all

        parameters:

        selection: list of int, slice object or list of bool
            Selection out of all the atoms that belong to the QM part.
        qmcalc: Calculator object
            QM-calculator.
        mmcalc1: Calculator object
            MM-calculator used for QM region.
        mmcalc2: Calculator object
            MM-calculator used for everything.
        vacuum: float or None
            Amount of vacuum to add around QM atoms.  Use None if QM
            calculator doesn't need a box.

        """
        self.selection = selection
        self.qmcalc = qmcalc
        self.mmcalc1 = mmcalc1
        self.mmcalc2 = mmcalc2
        self.vacuum = vacuum

        self.qmatoms = None
        self.center = None

        self.name = '{0}-{1}+{1}'.format(qmcalc.name, mmcalc1.name)

        Calculator.__init__(self)

    def initialize_qm(self, atoms):
        constraints = atoms.constraints
        atoms.constraints = []
        self.qmatoms = atoms[self.selection]
        atoms.constraints = constraints
        self.qmatoms.pbc = False
        if self.vacuum:
            self.qmatoms.center(vacuum=self.vacuum)
            self.center = self.qmatoms.positions.mean(axis=0)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.qmatoms is None:
            self.initialize_qm(atoms)

        self.qmatoms.positions = atoms.positions[self.selection]
        if self.vacuum:
            self.qmatoms.positions += (self.center -
                                       self.qmatoms.positions.mean(axis=0))

        energy = self.qmcalc.get_potential_energy(self.qmatoms)
        qmforces = self.qmcalc.get_forces(self.qmatoms)
        energy += self.mmcalc2.get_potential_energy(atoms)
        forces = self.mmcalc2.get_forces(atoms)

        if self.vacuum:
            qmforces -= qmforces.mean(axis=0)
        forces[self.selection] += qmforces

        energy -= self.mmcalc1.get_potential_energy(self.qmatoms)
        forces[self.selection] -= self.mmcalc1.get_forces(self.qmatoms)

        self.results['energy'] = energy
        self.results['forces'] = forces


class EIQMMM(Calculator):
    """Explicit interaction QMMM calculator."""
    implemented_properties = ['energy', 'forces']

    def __init__(self, selection, qmcalc, mmcalc, interaction,
                 vacuum=None, embedding=None, output=None):
        """EIQMMM object.

        The energy is calculated as::

                    _          _         _    _
            E = E  (R  ) + E  (R  ) + E (R  , R  )
                 QM  QM     MM  MM     I  QM   MM

        parameters:

        selection: list of int, slice object or list of bool
            Selection out of all the atoms that belong to the QM part.
        qmcalc: Calculator object
            QM-calculator.
        mmcalc: Calculator object
            MM-calculator.
        interaction: Interaction object
            Interaction between QM and MM regions.
        vacuum: float or None
            Amount of vacuum to add around QM atoms.  Use None if QM
            calculator doesn't need a box.
        embedding: Embedding object or None
            Specialized embedding object.  Use None in order to use the
            default one.
        output: None, '-', str or file-descriptor.
            File for logging information - default is no logging (None).

        """

        self.selection = selection

        self.qmcalc = qmcalc
        self.mmcalc = mmcalc
        self.interaction = interaction
        self.vacuum = vacuum
        self.embedding = embedding

        self.qmatoms = None
        self.mmatoms = None
        self.mask = None
        self.center = None  # center of QM atoms in QM-box

        self.name = '{0}+{1}+{2}'.format(qmcalc.name,
                                         interaction.name,
                                         mmcalc.name)

        self.output = convert_string_to_fd(output)

        Calculator.__init__(self)

    def initialize(self, atoms):
        self.mask = np.zeros(len(atoms), bool)
        self.mask[self.selection] = True

        constraints = atoms.constraints
        atoms.constraints = []  # avoid slicing of constraints
        self.qmatoms = atoms[self.mask]
        self.mmatoms = atoms[~self.mask]
        atoms.constraints = constraints

        self.qmatoms.pbc = False

        if self.vacuum:
            self.qmatoms.center(vacuum=self.vacuum)
            self.center = self.qmatoms.positions.mean(axis=0)
            print('Size of QM-cell after centering:',
                  self.qmatoms.cell.diagonal(), file=self.output)

        self.qmatoms.calc = self.qmcalc
        self.mmatoms.calc = self.mmcalc

        if self.embedding is None:
            self.embedding = Embedding()

        self.embedding.initialize(self.qmatoms, self.mmatoms)
        print('Embedding:', self.embedding, file=self.output)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.qmatoms is None:
            self.initialize(atoms)

        self.mmatoms.set_positions(atoms.positions[~self.mask])
        self.qmatoms.set_positions(atoms.positions[self.mask])

        if self.vacuum:
            shift = self.center - self.qmatoms.positions.mean(axis=0)
            self.qmatoms.positions += shift
        else:
            shift = (0, 0, 0)

        self.embedding.update(shift)

        ienergy, iqmforces, immforces = self.interaction.calculate(
            self.qmatoms, self.mmatoms, shift)

        qmenergy = self.qmatoms.get_potential_energy()
        mmenergy = self.mmatoms.get_potential_energy()
        energy = ienergy + qmenergy + mmenergy

        print('Energies: {0:12.3f} {1:+12.3f} {2:+12.3f} = {3:12.3f}'
              .format(ienergy, qmenergy, mmenergy, energy), file=self.output)

        qmforces = self.qmatoms.get_forces()
        mmforces = self.mmatoms.get_forces()

        mmforces += self.embedding.get_mm_forces()

        forces = np.empty((len(atoms), 3))
        forces[self.mask] = qmforces + iqmforces
        forces[~self.mask] = mmforces + immforces

        self.results['energy'] = energy
        self.results['forces'] = forces


def wrap(D, cell, pbc):
    """Wrap distances to nearest neighbor (minimum image convention)."""
    for i, periodic in enumerate(pbc):
        if periodic:
            d = D[:, i]
            L = cell[i]
            d[:] = (d + L / 2) % L - L / 2  # modify D inplace


class Embedding:
    def __init__(self, molecule_size=3, **parameters):
        """Point-charge embedding."""
        self.qmatoms = None
        self.mmatoms = None
        self.molecule_size = molecule_size
        self.virtual_molecule_size = None
        self.parameters = parameters

    def __repr__(self):
        return 'Embedding(molecule_size={0})'.format(self.molecule_size)

    def initialize(self, qmatoms, mmatoms):
        """Hook up embedding object to QM and MM atoms objects."""
        self.qmatoms = qmatoms
        self.mmatoms = mmatoms
        charges = mmatoms.calc.get_virtual_charges(mmatoms)
        self.pcpot = qmatoms.calc.embed(charges, **self.parameters)
        self.virtual_molecule_size = (self.molecule_size *
                                      len(charges) // len(mmatoms))

    def update(self, shift):
        """Update point-charge positions."""
        # Wrap point-charge positions to the MM-cell closest to the
        # center of the the QM box, but avoid ripping molecules apart:
        qmcenter = self.qmatoms.positions.mean(axis=0)
        # if counter ions are used, then molecule_size has more than 1 value
        if self.mmatoms.calc.name == 'combinemm':
            mask1 = self.mmatoms.calc.mask
            mask2 = ~mask1
            vmask1 = self.mmatoms.calc.virtual_mask
            vmask2 = ~vmask1
            apm1 = self.mmatoms.calc.apm1
            apm2 = self.mmatoms.calc.apm2
            spm1 = self.mmatoms.calc.atoms1.calc.sites_per_mol
            spm2 = self.mmatoms.calc.atoms2.calc.sites_per_mol
            pos = self.mmatoms.positions
            pos1 = pos[mask1].reshape((-1, apm1, 3))
            pos2 = pos[mask2].reshape((-1, apm2, 3))
            pos = (pos1, pos2)
        else:
            pos = (self.mmatoms.positions, )
            apm1 = self.molecule_size
            apm2 = self.molecule_size
            # This is only specific to calculators where apm != spm, 
            # i.e. TIP4P. Non-native MM calcs do not have this attr.
            if hasattr(self.mmatoms.calc, 'sites_per_mol'):
                spm1 = self.mmatoms.calc.sites_per_mol
                spm2 = self.mmatoms.calc.sites_per_mol
            else:
                spm1 = self.molecule_size
                spm2 = spm1
            mask1 = np.ones(len(self.mmatoms), dtype=bool)
            mask2 = mask1

        wrap_pos = np.zeros_like(self.mmatoms.positions)
        com_all = []
        apm = (apm1, apm2)
        mask = (mask1, mask2)
        spm = (spm1, spm2)
        for p, n, m, vn in zip(pos, apm, mask, spm):
            positions = p.reshape((-1, n, 3)) + shift

            # Distances from the center of the QM box to the first atom of
            # each molecule:
            distances = positions[:, 0] - qmcenter

            wrap(distances, self.mmatoms.cell.diagonal(), self.mmatoms.pbc)
            offsets = distances - positions[:, 0]
            positions += offsets[:, np.newaxis] + qmcenter

            # Geometric center positions for each mm mol for LR cut
            com = np.array([p.mean(axis=0) for p in positions])
            # Need per atom for C-code:
            com_pv = np.repeat(com, vn, axis=0)
            com_all.append(com_pv)

            wrap_pos[m] = positions.reshape((-1, 3))

        positions = wrap_pos.copy()
        positions = self.mmatoms.calc.add_virtual_sites(positions)

        if self.mmatoms.calc.name == 'combinemm':
            com_pv = np.zeros_like(positions)
            for ii, m in enumerate((vmask1, vmask2)):
                com_pv[m] = com_all[ii]

        # compatibility with gpaw versions w/o LR cut in PointChargePotential
        if 'rc2' in self.parameters:
            self.pcpot.set_positions(positions, com_pv=com_pv)
        else:
            self.pcpot.set_positions(positions)

    def get_mm_forces(self):
        """Calculate the forces on the MM-atoms from the QM-part."""
        f = self.pcpot.get_forces(self.qmatoms.calc)
        return self.mmatoms.calc.redistribute_forces(f)


def combine_lj_lorenz_berthelot(sigmaqm, sigmamm,
                                epsilonqm, epsilonmm):
    """Combine LJ parameters according to the Lorenz-Berthelot rule"""
    sigma = []
    epsilon = []
    # check if input is tuple of vals for more than 1 mm calc, or only for 1.
    if type(sigmamm) == tuple:
        numcalcs = len(sigmamm)
    else:
        numcalcs = 1  # if only 1 mm calc, eps and sig are simply np arrays
        sigmamm = (sigmamm, )
        epsilonmm = (epsilonmm, )
    for cc in range(numcalcs):
        sigma_c = np.zeros((len(sigmaqm), len(sigmamm[cc])))
        epsilon_c = np.zeros_like(sigma_c)

        for ii in range(len(sigmaqm)):
            sigma_c[ii, :] = (sigmaqm[ii] + sigmamm[cc]) / 2
            epsilon_c[ii, :] = (epsilonqm[ii] * epsilonmm[cc])**0.5
        sigma.append(sigma_c)
        epsilon.append(epsilon_c)

    if numcalcs == 1:  # retain original, 1 calc function
        sigma = np.array(sigma[0])
        epsilon = np.array(epsilon[0])

    return sigma, epsilon


class LJInteractionsGeneral:
    name = 'LJ-general'

    def __init__(self, sigmaqm, epsilonqm, sigmamm, epsilonmm,
                 qm_molecule_size, mm_molecule_size=3,
                 rc=np.Inf, width=1.0):
        """General Lennard-Jones type explicit interaction.

        sigmaqm: array
            Array of sigma-parameters which should have the length of the QM
            subsystem
        epsilonqm: array
            As sigmaqm, but for epsilon-paramaters
        sigmamm: Either array (A) or tuple (B)
            A (no counterions): 
                Array of sigma-parameters with the length of the smallest 
                repeating atoms-group (i.e. molecule) of the MM subsystem
            B (counterions):
                Tuple: (arr1, arr2), where arr1 is an array of sigmas with
                the length of counterions in the MM subsystem, and 
                arr2 is the array from A.
        epsilonmm: array or tuple
            As sigmamm but for epsilon-parameters.
        qm_molecule_size: int
            number of atoms of the smallest repeating atoms-group (i.e. 
            molecule) in the QM subystem (often just the number of atoms 
            in the QM subsystem)
        mm_molecule_size: int
            as qm_molecule_size but for the MM subsystem. Will be overwritten
            if counterions are present in the MM subsystem (via the CombineMM
            calculator)

        """
        self.sigmaqm = sigmaqm
        self.epsilonqm = epsilonqm
        self.sigmamm = sigmamm
        self.epsilonmm = epsilonmm
        self.qms = qm_molecule_size
        self.mms = mm_molecule_size
        self.rc = rc
        self.width = width
        self.combine_lj()

    def combine_lj(self):
        self.sigma, self.epsilon = combine_lj_lorenz_berthelot(
            self.sigmaqm, self.sigmamm, self.epsilonqm, self.epsilonmm)

    def calculate(self, qmatoms, mmatoms, shift):
        epsilon = self.epsilon
        sigma = self.sigma

        # loop over possible multiple mm calculators
        # currently 1 or 2, but could be generalized in the future...
        apm1 = self.mms
        mask1 = np.ones(len(mmatoms), dtype=bool)
        mask2 = mask1
        apm = (apm1, )
        sigma = (sigma, )
        epsilon = (epsilon, )
        if hasattr(mmatoms.calc, 'name'):
            if mmatoms.calc.name == 'combinemm':
                mask1 = mmatoms.calc.mask
                mask2 = ~mask1
                apm1 = mmatoms.calc.apm1
                apm2 = mmatoms.calc.apm2
                apm = (apm1, apm2)
                sigma = sigma[0]  # Was already loopable 2-tuple
                epsilon = epsilon[0]

        mask = (mask1, mask2)
        e_all = 0
        qmforces_all = np.zeros_like(qmatoms.positions)
        mmforces_all = np.zeros_like(mmatoms.positions)

        # zip stops at shortest tuple so we dont double count
        # case_koopmanss of no counter ions.
        for n, m, eps, sig in zip(apm, mask, epsilon, sigma):
            mmpositions = self.update(qmatoms, mmatoms[m], n, shift)
            qmforces = np.zeros_like(qmatoms.positions)
            mmforces = np.zeros_like(mmatoms[m].positions)
            energy = 0.0

            qmpositions = qmatoms.positions.reshape((-1, self.qms, 3))

            for q, qmpos in enumerate(qmpositions):  # molwise loop
                # cutoff from first atom of each mol
                R00 = mmpositions[:, 0] - qmpos[0, :]
                d002 = (R00**2).sum(1)
                d00 = d002**0.5
                x1 = d00 > self.rc - self.width
                x2 = d00 < self.rc
                x12 = np.logical_and(x1, x2)
                y = (d00[x12] - self.rc + self.width) / self.width
                t = np.zeros(len(d00))
                t[x2] = 1.0
                t[x12] -= y**2 * (3.0 - 2.0 * y)
                dt = np.zeros(len(d00))
                dt[x12] -= 6.0 / self.width * y * (1.0 - y)
                for qa in range(len(qmpos)):
                    if ~np.any(eps[qa, :]):
                        continue
                    R = mmpositions - qmpos[qa, :]
                    d2 = (R**2).sum(2)
                    c6 = (sig[qa, :]**2 / d2)**3
                    c12 = c6**2
                    e = 4 * eps[qa, :] * (c12 - c6)
                    energy += np.dot(e.sum(1), t)
                    f = t[:, None, None] * (24 * eps[qa, :] *
                         (2 * c12 - c6) / d2)[:, :, None] * R
                    f00 = - (e.sum(1) * dt / d00)[:, None] * R00
                    mmforces += f.reshape((-1, 3))
                    qmforces[q * self.qms + qa, :] -= f.sum(0).sum(0)
                    qmforces[q * self.qms, :] -= f00.sum(0)
                    mmforces[::n, :] += f00

                e_all += energy
                qmforces_all += qmforces
                mmforces_all[m] += mmforces

        return e_all, qmforces_all, mmforces_all

    def update(self, qmatoms, mmatoms, n, shift):
        # Wrap point-charge positions to the MM-cell closest to the
        # center of the the QM box, but avoid ripping molecules apart:
        qmcenter = qmatoms.cell.diagonal() / 2
        positions = mmatoms.positions.reshape((-1, n, 3)) + shift

        # Distances from the center of the QM box to the first atom of
        # each molecule:
        distances = positions[:, 0] - qmcenter

        wrap(distances, mmatoms.cell.diagonal(), mmatoms.pbc)
        offsets = distances - positions[:, 0]
        positions += offsets[:, np.newaxis] + qmcenter

        return positions


class LJInteractions:
    name = 'LJ'

    def __init__(self, parameters):
        """Lennard-Jones type explicit interaction.

        parameters: dict
            Mapping from pair of atoms to tuple containing epsilon and sigma
            for that pair.

        Example:

            lj = LJInteractions({('O', 'O'): (eps, sigma)})

        """
        self.parameters = {}
        for (symbol1, symbol2), (epsilon, sigma) in parameters.items():
            Z1 = atomic_numbers[symbol1]
            Z2 = atomic_numbers[symbol2]
            self.parameters[(Z1, Z2)] = epsilon, sigma
            self.parameters[(Z2, Z1)] = epsilon, sigma

    def calculate(self, qmatoms, mmatoms, shift):
        qmforces = np.zeros_like(qmatoms.positions)
        mmforces = np.zeros_like(mmatoms.positions)
        species = set(mmatoms.numbers)
        energy = 0.0
        for R1, Z1, F1 in zip(qmatoms.positions, qmatoms.numbers, qmforces):
            for Z2 in species:
                if (Z1, Z2) not in self.parameters:
                    continue
                epsilon, sigma = self.parameters[(Z1, Z2)]
                mask = (mmatoms.numbers == Z2)
                D = mmatoms.positions[mask] + shift - R1
                wrap(D, mmatoms.cell.diagonal(), mmatoms.pbc)
                d2 = (D**2).sum(1)
                c6 = (sigma**2 / d2)**3
                c12 = c6**2
                energy += 4 * epsilon * (c12 - c6).sum()
                f = 24 * epsilon * ((2 * c12 - c6) / d2)[:, np.newaxis] * D
                F1 -= f.sum(0)
                mmforces[mask] += f
        return energy, qmforces, mmforces


class RescaledCalculator(Calculator):
    """Rescales length and energy of a calculators to match given
    lattice constant and bulk modulus

    Useful for MM calculator used within a :class:`ForceQMMM` model.
    See T. D. Swinburne and J. R. Kermode, Phys. Rev. B 96, 144102 (2017)
    for a derivation of the scaling constants.
    """
    implemented_properties = ['forces', 'energy', 'stress']

    def __init__(self, mm_calc,
                 qm_lattice_constant, qm_bulk_modulus,
                 mm_lattice_constant, mm_bulk_modulus):
        Calculator.__init__(self)
        self.mm_calc = mm_calc
        self.alpha = qm_lattice_constant / mm_lattice_constant
        self.beta = mm_bulk_modulus / qm_bulk_modulus / (self.alpha**3)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        # mm_pos = atoms.get_positions()
        scaled_atoms = atoms.copy()

        # scaled_atoms.positions = mm_pos/self.alpha
        mm_cell = atoms.get_cell()
        scaled_atoms.set_cell(mm_cell / self.alpha, scale_atoms=True)

        forces = self.mm_calc.get_forces(scaled_atoms)
        energy = self.mm_calc.get_potential_energy(scaled_atoms)
        stress = self.mm_calc.get_stress(scaled_atoms)

        self.results = {'energy': energy / self.beta,
                        'forces': forces / (self.beta * self.alpha),
                        'stress': stress / (self.beta * self.alpha**3)}


class ForceConstantCalculator(Calculator):
    """
    Compute forces base_koopmansd on provided force-constant matrix

    Useful with `ForceQMMM` to do harmonic QM/MM using force constants
    of QM method.
    """
    implemented_properties = ['forces', 'energy']

    def __init__(self, D, ref, f0):
        """
        Parameters:

        D: matrix or sparse matrix, shape `(3*len(ref), 3*len(ref))`
            Force constant matrix.
            Sign convention is `D_ij = d^2E/(dx_i dx_j), so
            `force = -D.dot(displacement)`
        ref: ase_koopmans.atoms.Atoms
            Atoms object for reference configuration
        f0: array, shape `(len(ref), 3)`
            Value of forces at reference configuration
        """
        assert D.shape[0] == D.shape[1]
        assert D.shape[0] // 3 == len(ref)
        self.D = D
        self.ref = ref
        self.f0 = f0
        self.size = len(ref)
        Calculator.__init__(self)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        u = atoms.positions - self.ref.positions
        f = -self.D.dot(u.reshape(3 * self.size))
        forces = np.zeros((len(atoms), 3))
        forces[:, :] = f.reshape(self.size, 3)
        self.results['forces'] = forces + self.f0
        self.results['energy'] = 0.0


class ForceQMMM(Calculator):
    """
    Force-base_koopmansd QM/MM calculator

    QM forces are computed using a buffer region and then mixed abruptly
    with MM forces:

        F^i_QMMM = {   F^i_QM    if i in QM region
                   {   F^i_MM    otherwise

    cf. N. Bernstein, J. R. Kermode, and G. Csanyi,
    Rep. Prog. Phys. 72, 026501 (2009)
    and T. D. Swinburne and J. R. Kermode, Phys. Rev. B 96, 144102 (2017).
    """
    implemented_properties = ['forces', 'energy']

    def __init__(self,
                 atoms,
                 qm_selection_mask,
                 qm_calc,
                 mm_calc,
                 buffer_width,
                 vacuum=5.,
                 zero_mean=True):
        """
        ForceQMMM calculator

        Parameters:

        qm_selection_mask: list of ints, slice object or bool list/array
            Selection out of atoms that belong to the QM region.
        qm_calc: Calculator object
            QM-calculator.
        mm_calc: Calculator object
            MM-calculator (should be scaled, see :class:`RescaledCalculator`)
            Can use `ForceConstantCalculator` base_koopmansd on QM force constants, if
            available.
        vacuum: float or None
            Amount of vacuum to add around QM atoms.
        zero_mean: bool
            If True, add a correction to zero the mean force in each direction
        """

        if len(atoms[qm_selection_mask]) == 0:
            raise ValueError("no QM atoms selected!")

        self.qm_selection_mask = qm_selection_mask
        self.qm_calc = qm_calc
        self.mm_calc = mm_calc
        self.vacuum = vacuum
        self.buffer_width = buffer_width
        self.zero_mean = zero_mean

        self.qm_buffer_mask = None
        self.cell = None
        self.qm_shift = None

        Calculator.__init__(self)

    def initialize_qm_buffer_mask(self, atoms):
        """
        Initialises system to perform qm calculation
        """

        # get the radius of the qm_selection in non periodic directions
        qm_positions = atoms[self.qm_selection_mask].get_positions()
        # identify qm radius as an larges distance from the center
        # of the cluster (overestimation)
        qm_center = qm_positions.mean(axis=0)

        non_pbc_directions = np.logical_not(self.atoms.pbc)

        centered_positions = atoms.get_positions()

        for i, non_pbc in enumerate(non_pbc_directions):
            if non_pbc:
                qm_positions.T[i] -= qm_center[i]
                centered_positions.T[i] -= qm_center[i]

        qm_radius = np.linalg.norm(qm_positions.T, axis=1).max()
        self.cell = self.atoms.cell.copy()

        for i, non_pbc in enumerate(non_pbc_directions):
            if non_pbc:
                self.cell[i][i] = 2.0 * (qm_radius +
                                         self.buffer_width +
                                         self.vacuum)

        # identify atoms in region < qm_radius + buffer
        distances_from_center = np.linalg.norm(
            centered_positions.T[non_pbc_directions].T, axis=1)

        self.qm_buffer_mask = (distances_from_center <
                               qm_radius + self.buffer_width)

        # exclude atoms that are too far (in case_koopmans of non spherical region)
        for i, buffer_atom in enumerate(self.qm_buffer_mask &
                                        np.logical_not(self.qm_selection_mask)):
            if buffer_atom:
                distance = np.linalg.norm(
                    (qm_positions -
                     centered_positions[i]).T[non_pbc_directions].T, axis=1)
                if distance.min() > self.buffer_width:
                    self.qm_buffer_mask[i] = False

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        if self.qm_buffer_mask is None:
            self.initialize_qm_buffer_mask(atoms)

        # initialize the object
        # qm_buffer_atoms = atoms.copy()
        qm_buffer_atoms = atoms[self.qm_buffer_mask]
        del qm_buffer_atoms.constraints

        qm_buffer_atoms.set_cell(self.cell)
        qm_shift = (0.5 * qm_buffer_atoms.cell.diagonal() -
                    qm_buffer_atoms.positions.mean(axis=0))

        qm_buffer_atoms.set_cell(self.cell)
        qm_buffer_atoms.positions += qm_shift

        forces = self.mm_calc.get_forces(atoms)

        qm_forces = self.qm_calc.get_forces(qm_buffer_atoms)
        forces[self.qm_selection_mask] = \
            qm_forces[self.qm_selection_mask[self.qm_buffer_mask]]

        if self.zero_mean:
            # Target is that: forces.sum(axis=1) == [0., 0., 0.]
            forces[:] -= forces.mean(axis=0)

        self.results['forces'] = forces
        self.results['energy'] = 0.0

