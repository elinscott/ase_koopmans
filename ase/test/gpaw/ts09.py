import unittest
import numpy as np

import ase.io
from ase.calculators.vdwcorrection import vdWTkatchenko09prl
from ase.parallel import barrier, parprint
from ase.build import molecule

try:
    from gpaw import GPAW
except ImportError:
    # Skip test if GPAW installation is broken:
    raise unittest.SkipTest

from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw.analyse.vdwradii import vdWradii
from gpaw.cluster import Cluster
from gpaw.test import equal

h = 0.4
s = Cluster(molecule('LiH'))
s.minimal_box(3., h=h)


def print_charge_and_check(hp, q=0, label='unpolarized'):
    q_a = np.array(hp.get_charges())
    parprint('Charges ({0})='.format(label), q_a, ', sum=', q_a.sum())
    equal(q_a.sum(), q, 0.03)
    return q_a

# spin unpolarized

if 1:
    out_traj = 'LiH.traj'
    out_txt = 'LiH.txt'

    cc = GPAW(h=h, xc='PBE', txt=out_txt)
    hp = HirshfeldPartitioning(cc)
    c = vdWTkatchenko09prl(hp,
                           vdWradii(s.get_chemical_symbols(), 'PBE'))
    s.set_calculator(c)
    E = s.get_potential_energy()
    F_ac = s.get_forces()
    s.write(out_traj)
    q_a = print_charge_and_check(hp)

    barrier()

    # test I/O, accuracy due to text output
    accuracy = 1.e-5
    for fname in [out_traj, out_txt]:
        s_out = ase.io.read(fname)
        equal(s_out.get_potential_energy(), E, accuracy)
        for fi, fo in zip(F_ac, s_out.get_forces()):
            equal(fi, fo, accuracy)
        if fname.endswith('.traj'):
            assert(s_out.get_calculator().parameters['calculator'] == 'gpaw')


# spin polarized

if 0:
    ccs = GPAW(h=h, xc='PBE', spinpol=True,
               txt=None)
    hps = HirshfeldPartitioning(ccs)
    cs = vdWTkatchenko09prl(hps, vdWradii(s.get_chemical_symbols(), 'PBE'))
    s.set_calculator(cs)
    Es = s.get_potential_energy()
    Fs_ac = s.get_forces()

    qs_a = print_charge_and_check(hps, label='spin')

    equal(q_a, qs_a, 1.e-6)
    equal(E, Es, 1.e-4)
    equal(F_ac, Fs_ac, 1.e-4)

# charged

if 0:
    cc.set(charge=1)
    hpp = HirshfeldPartitioning(cc)
    cp = vdWTkatchenko09prl(hpp,
                            vdWradii(s.get_chemical_symbols(), 'PBE'))
    s.set_calculator(cp)
    E = s.get_potential_energy()
    
    print_charge_and_check(hpp, 1, label='+1')
