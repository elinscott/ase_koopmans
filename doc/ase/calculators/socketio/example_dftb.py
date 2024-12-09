import sys
from ase_koopmans.build import molecule
from ase_koopmans.calculators.dftb import Dftb
from ase_koopmans.calculators.socketio import SocketIOCalculator
from ase_koopmans.optimize import BFGS


atoms = molecule('H2O')
dftb = Dftb(Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='"p"',
            Hamiltonian_MaxAngularMomentum_H='"s"',
            Driver_='',
            Driver_Socket_='',
            Driver_Socket_File='Hello')
opt = BFGS(atoms, trajectory='test.traj')


with SocketIOCalculator(dftb, log=sys.stdout, unixsocket='Hello') as calc:
    atoms.calc = calc
    opt.run(fmax=0.01)
