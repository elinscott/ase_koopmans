"""Langevin dynamics class."""

import numpy as np

from ase.md.md import MolecularDynamics
from ase.parallel import world


class Langevin(MolecularDynamics):
    """Langevin (constant N, V, T) molecular dynamics.

    Usage: Langevin(atoms, dt, temperature, friction)

    atoms
        The list of atoms.

    dt
        The time step.

    temperature
        The desired temperature, in energy units.

    friction
        A friction coefficient, typically 1e-4 to 1e-2.

    selectlinear
        Selection of atoms that are not part of rigid linear triatomic
        molecules. Default is None, which means that the atoms object
        should not include molecules kept rigid using FixLinearTriatomic
        constraints. If these are present, "selectlinear" must be specified
        and a sequence ABCABC... of atoms of linear triatomic molecules of
        the same type is assumed.

    fixcm
        If True, the position and momentum of the center of mass is
        kept unperturbed.  Default: True.

    rng
        Random number generator, by default numpy.random.  Must have a
        standard_normal method matching the signature of
        numpy.random.standard_normal.

    The temperature and friction are normally scalars, but in principle one
    quantity per atom could be specified by giving an array.

    RATTLE constraints can be used with these propagators, see:
    E. V.-Eijnden, and G. Ciccotti, Chem. Phys. Lett. 429, 310 (2006)

    The propagator is Equation 23 (Eq. 39 if RATTLE constraints are used)
    of the above reference.  That reference also contains another
    propagator in Eq. 21/34; but that propagator is not quasi-symplectic
    and gives a systematic offset in the temperature at large time steps.

    This dynamics accesses the atoms using Cartesian coordinates."""

    # Helps Asap doing the right thing.  Increment when changing stuff:
    _lgv_version = 3

    def __init__(self, atoms, timestep, temperature, friction,
                 selectlinear=None, fixcm=True, trajectory=None, logfile=None,
                 loginterval=1, communicator=world, rng=np.random):
        for constraint in atoms.constraints:
            if (constraint.todict()['name'] == 'FixLinearTriatomic' and
               selectlinear is None):
                raise ValueError('Specify "selectlinear" when using'
                                 'FixLinearTriatomic constraints')
        self.temp = temperature
        self.fr = friction
        self.selectlinear = selectlinear
        self.fixcm = fixcm  # will the center of mass be held fixed?
        self.communicator = communicator
        self.rng = rng

        MolecularDynamics.__init__(self, atoms, timestep, trajectory,
                                   logfile, loginterval)

        if self.selectlinear is not None:
            self.mask = np.zeros(len(atoms), bool)
            self.mask[self.selectlinear] = True

            # Friction on central atom must be zero
            if np.isscalar(self.fr):
                self.fr = np.full((len(atoms), 3), self.fr)
                frl = self.fr[~self.mask]
                frl[1::3] = 0.0
                self.fr[~self.mask] = frl
            else:
                frl = self.fr[~self.mask]
                frl[1::3] = 0.0
                self.fr[~self.mask] = frl

            masseslinear = self.masses[~self.mask]
            m_a = masseslinear[0]
            m_b = masseslinear[1]
            m_c = masseslinear[2]
            atomslinear = atoms[~self.mask]
            r_ab = atomslinear.get_distance(0, 1)
            r_bc = atomslinear.get_distance(1, 2)
            r_ac = r_ab + r_bc
            self.m_ab = m_a * m_b
            self.m_bc = m_b * m_c
            self.m_ac = m_a * m_c
            self.c_a = r_bc / r_ac
            self.c_c = r_ab / r_ac
            self.mr_bc = (m_b / m_c)**0.5
            self.mr_ba = (m_b / m_a)**0.5
            self.mr_ca = (m_c / m_a)**0.5
            self.mr_ac = (m_a / m_c)**0.5
            self.n_a = self.c_a / (self.c_a**2 * self.m_bc +
                                   self.c_c**2 * self.m_ab + self.m_ac)
            self.n_c = self.c_c / (self.c_a**2 * self.m_bc +
                                   self.c_c**2 * self.m_ab + self.m_ac)

        self.updatevars()

    def todict(self):
        d = MolecularDynamics.todict(self)
        d.update({'temperature': self.temp,
                  'friction': self.fr,
                  'fix-cm': self.fixcm})
        return d

    def set_temperature(self, temperature):
        self.temp = temperature
        self.updatevars()

    def set_friction(self, friction):
        self.fr = friction
        self.updatevars()

    def set_timestep(self, timestep):
        self.dt = timestep
        self.updatevars()

    def updatevars(self):
        dt = self.dt
        T = self.temp
        fr = self.fr
        print(fr)
        masses = self.masses
        sigma = np.sqrt(2 * T * fr / masses)

        self.c1 = dt / 2. - dt * dt * fr / 8.
        self.c2 = dt * fr / 2 - dt * dt * fr * fr / 8.
        self.c3 = np.sqrt(dt) * sigma / 2. - dt**1.5 * fr * sigma / 8.
        self.c5 = dt**1.5 * sigma / (2 * np.sqrt(3))
        self.c4 = fr / 2. * self.c5

        # Works in parallel Asap, #GLOBAL number of atoms:
        self.natoms = self.atoms.get_number_of_atoms()

    def step(self, f):
        atoms = self.atoms
        natoms = len(atoms)

        # This velocity as well as xi, eta and a few other variables are stored
        # as attributes, so Asap can do its magic when atoms migrate between
        # processors.
        self.v = atoms.get_velocities()

        self.xi = self.rng.standard_normal(size=(natoms, 3))
        self.eta = self.rng.standard_normal(size=(natoms, 3))

        if self.selectlinear is not None:
            self.xi, self.eta = self.redistribute()

        if self.communicator is not None:
            self.communicator.broadcast(self.xi, 0)
            self.communicator.broadcast(self.eta, 0)

        # First halfstep in the velocity.
        self.v += (self.c1 * f / self.masses - self.c2 * self.v +
                   self.c3 * self.xi - self.c4 * self.eta)

        # Full step in positions
        x = atoms.get_positions()
        if self.fixcm:
            old_cm = atoms.get_center_of_mass()
        # Step: x^n -> x^(n+1) - this applies constraints if any.
        atoms.set_positions(x + self.dt * self.v + self.c5 * self.eta)
        if self.fixcm:
            new_cm = atoms.get_center_of_mass()
            d = old_cm - new_cm
            # atoms.translate(d)  # Does not respect constraints
            atoms.set_positions(atoms.get_positions() + d)

        # recalc velocities after RATTLE constraints are applied
        self.v = (self.atoms.get_positions() - x -
                  self.c5 * self.eta) / self.dt
        f = atoms.get_forces(md=True)

        # Update the velocities
        self.v += (self.c1 * f / self.masses - self.c2 * self.v +
                   self.c3 * self.xi - self.c4 * self.eta)

        if self.fixcm:  # subtract center of mass vel
            v_cm = self._get_com_velocity()
            self.v -= v_cm

        # Second part of RATTLE taken care of here
        atoms.set_momenta(self.v * self.masses)

        return f

    def redistribute(self):
        m_ab = self.m_ab
        m_bc = self.m_bc
        m_ac = self.m_ac
        c_a = self.c_a
        c_c = self.c_c
        mr_bc = self.mr_bc
        mr_ba = self.mr_ba
        mr_ca = self.mr_ca
        mr_ac = self.mr_ac
        n_a = self.n_a
        n_c = self.n_c

        xi = self.xi
        eta = self.eta
        xilin = xi[~self.mask]
        etalin = eta[~self.mask]
        xinew = np.zeros_like(xi)
        etanew = np.zeros_like(eta)
        xir = np.zeros_like(xilin)
        etar = np.zeros_like(etalin)

        # Avoid redistributing for atoms in selectlinear
        xinew[self.mask] = xi[self.mask]
        etanew[self.mask] = eta[self.mask]

        xir[::3, :] = ((1 - n_a * m_bc * c_a) * xilin[::3, :] -
                       n_a * (m_ab * c_c * mr_ca * xilin[2::3, :] -
                       m_ac * mr_ba * xilin[1::3, :]))
        etar[::3, :] = ((1 - n_a * m_bc * c_a) * etalin[::3, :] -
                        n_a * (m_ab * c_c * mr_ca * etalin[2::3, :] -
                        m_ac * mr_ba * etalin[1::3, :]))
        xir[2::3, :] = ((1 - n_c * m_ab * c_c) * xilin[2::3, :] -
                        n_c * (m_bc * c_a * mr_ac * xilin[::3, :] -
                        m_ac * mr_bc * xilin[1::3, :]))
        etar[2::3, :] = ((1 - n_c * m_ab * c_c) * etalin[2::3, :] -
                         n_c * (m_bc * c_a * mr_ac * etalin[::3, :] -
                         m_ac * mr_bc * etalin[1::3, :]))

        xinew[~self.mask] = xir
        etanew[~self.mask] = etar

        return xinew, etanew

    def _get_com_velocity(self):
        """Return the center of mass velocity.

        Internal use only.  This function can be reimplemented by Asap.
        """
        return np.dot(self.masses.flatten(), self.v) / self.masses.sum()
