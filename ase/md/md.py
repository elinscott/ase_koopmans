"""Molecular Dynamics."""

import warnings
import numpy as np

from ase.optimize.optimize import Dynamics
from ase.md.logger import MDLogger


class MolecularDynamics(Dynamics):
    """Base-class for all MD classes."""
    def __init__(self, atoms, timestep, trajectory, logfile=None,
                 loginterval=1, append_trajectory=False):

        # dt as to be attached _before_ parent class is initialized
        self.dt = timestep

        Dynamics.__init__(self, atoms, logfile=None, trajectory=trajectory,
                          append_trajectory=append_trajectory)

        self.masses = self.atoms.get_masses()
        self.max_steps = None

        if 0 in self.masses:
            warnings.warn('Zero mass encountered in atoms; this will '
                          'likely lead to errors if the massless atoms '
                          'are unconstrained.')

        self.masses.shape = (-1, 1)

        if not self.atoms.has('momenta'):
            self.atoms.set_momenta(np.zeros([len(self.atoms), 3]))

        if logfile:
            self.attach(MDLogger(dyn=self, atoms=atoms, logfile=logfile),
                        interval=loginterval)

    def todict(self):
        return {'type': 'molecular-dynamics',
                'md-type': self.__class__.__name__,
                'timestep': self.dt}

    def irun(self, steps=50):
        """ Call Dynamics.irun and adjust max_steps """
        self.max_steps = steps + self.nsteps
        return Dynamics.irun(self)

    def run(self, steps=50):
        """ Call Dynamics.run and adjust max_steps """
        self.max_steps = steps + self.nsteps
        return Dynamics.run(self)

    def get_time(self):
        return self.nsteps * self.dt

    def converged(self):
        """ MD is 'converged' when number of maximum steps is reached. """
        return self.nsteps >= self.max_steps
