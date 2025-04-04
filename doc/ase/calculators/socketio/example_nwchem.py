import sys

from ase_koopmans.build import molecule
from ase_koopmans.calculators.nwchem import NWChem
from ase_koopmans.calculators.socketio import SocketIOCalculator
from ase_koopmans.optimize import BFGS

atoms = molecule('H2O')
atoms.rattle(stdev=0.1)

unixsocket = 'ase_nwchem'

nwchem = NWChem(theory='scf',
                task='optimize',
                driver={'socket': {'unix': unixsocket}})

opt = BFGS(atoms, trajectory='opt.traj',
           logfile='opt.log')

with SocketIOCalculator(nwchem, log=sys.stdout,
                        unixsocket=unixsocket) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
