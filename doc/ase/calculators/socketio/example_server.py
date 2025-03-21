import sys

from ase_koopmans.build import molecule
from ase_koopmans.calculators.socketio import SocketIOCalculator
from ase_koopmans.io import write
from ase_koopmans.optimize import BFGS

unixsocket = 'ase_server_socket'

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)
write('initial.traj', atoms)

opt = BFGS(atoms, trajectory='opt.driver.traj', logfile='opt.driver.log')

with SocketIOCalculator(log=sys.stdout,
                        unixsocket=unixsocket) as calc:
    # Server is now running and waiting for connections.
    # If you want to launch the client process here directly,
    # instead of manually in the terminal, uncomment these lines:
    #
    # from subprocess import Popen
    # proc = Popen([sys.executable, 'example_client_gpaw.py'])

    atoms.calc = calc
    opt.run(fmax=0.05)
