from ase import Atoms, Atom
from ase.vibrations.placzek import Placzek, Profeta
from ase.vibrations.albrecht import Albrecht

from gpaw import GPAW
from gpaw.lrtddft.kssingle import KSSingles
from gpaw.test import equal
from gpaw.analyse.overlap import Overlap

"""
Test Placzek Raman implementation
"""

txt = '-'
txt = None
load = True
load = False
xc = 'LDA'

# ---------------------------------------------------
# run

R = 0.7  # approx. experimental bond length
a = 4.0
c = 5.0
H2 = Atoms([Atom('H', (a / 2, a / 2, (c - R) / 2)),
            Atom('H', (a / 2, a / 2, (c + R) / 2))],
           cell=(a, a, c))

gsname = exname = 'rraman'
exkwargs={'eps':0.0, 'jend':3}

if 1:
    calc = GPAW(xc=xc, nbands=7,
                convergence={'bands':3},
                spinpol=False, #eigensolver='rmm-diis',
                txt=txt)
    H2.set_calculator(calc)
    #H2.get_potential_energy()

    pz = Placzek(H2, KSSingles, gsname=gsname, exname=exname,
                 exkwargs=exkwargs,
                 # XXX full does not work in parallel due to boxes
                 # on different nodes
                 #overlap=lambda x, y: Overlap(x).full(y)[0],
                 overlap=lambda x, y: Overlap(x).pseudo(y)[0],
                 txt=txt)
    pz.run()

# ---------------------------------------------------
# check

"""Different Placzeck implementations should agree"""

om = 5
pz = Placzek(H2, KSSingles, gsname=gsname, exname=exname, txt=txt)
pzi = pz.absolute_intensity(omega=om)[-1]

pr = Profeta(H2, KSSingles, gsname=gsname, exname=exname,
             approximation='Placzek', txt=txt)
pri = pr.absolute_intensity(omega=om)[-1]
equal(pzi, pri, 0.1)

pr = Profeta(H2, KSSingles, gsname=gsname, exname=exname, 
             overlap=True,
             approximation='Placzek', txt=txt)
pri = pr.absolute_intensity(omega=om)[-1]
equal(pzi, pri, 0.1)

"""Albrecht and Placzek are approximately equal"""

al = Albrecht(H2, KSSingles, gsname=gsname, exname=exname, 
             overlap=True,
             approximation='Albrecht', txt=txt)
ali = al.absolute_intensity(omega=om)[-1]
equal(pzi, ali, 1.5)

"""Albrecht A and P-P are approximately equal"""

pr = Profeta(H2, KSSingles, gsname=gsname, exname=exname, 
             overlap=True,
             approximation='P-P', txt=txt)
pri = pr.absolute_intensity(omega=om)[-1]

al = Albrecht(H2, KSSingles, gsname=gsname, exname=exname, 
             overlap=True,
             approximation='Albrecht A', txt=txt)
ali = al.absolute_intensity(omega=om)[-1]
equal(pri, ali, 3)

"""Albrecht B+C and Profeta are approximately equal"""

pr = Profeta(H2, KSSingles, gsname=gsname, exname=exname, 
             overlap=True,
             approximation='Profeta', txt=txt)
pri = pr.absolute_intensity(omega=om)[-1]

al = Albrecht(H2, KSSingles, gsname=gsname, exname=exname, 
             overlap=True,
             approximation='Albrecht BC', txt=txt)
ali = al.absolute_intensity(omega=om)[-1]
equal(pri, ali, 3)

