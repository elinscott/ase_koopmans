from __future__ import print_function
from ase_koopmans.io import read
from ase_koopmans.calculators.socketio import SocketClient
from gpaw import GPAW, Mixer

# The atomic numbers are not transferred over the socket, so we have to
# read the file
atoms = read('initial.traj')
unixsocket = 'ase_server_socket'

atoms.calc = GPAW(mode='lcao',
                  basis='dzp',
                  txt='gpaw.client.txt',
                  mixer=Mixer(0.7, 7, 20.0))

client = SocketClient(unixsocket=unixsocket)

# Each step of the loop changes the atomic positions, but the generator
# yields None.
for i, _ in enumerate(client.irun(atoms, use_stress=False)):
    print('step:', i)
