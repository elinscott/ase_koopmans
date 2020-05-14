import numpy as np

from ase import Atoms, Atom

ase = """\
 H   HH HHH
H H H   H
HHH  H  HH
H H   H H
H H HH  HHH"""

d = 1.2

logo = Atoms()
for i, line in enumerate(ase.split('\n')):
    for j, c in enumerate(line):
        if c == 'H':
            logo.append(Atom('H', [d * j, d * i, 0]))
logo.set_cell((15, 15, 2))
logo.center()

if 1:
    from gpaw import GPAW
    calc = GPAW()
    logo.calc = calc
    e = logo.get_potential_energy()
    calc.write('logo2.gpw')
if 0:
    from gpaw import GPAW
    calc = GPAW('logo2.gpw', idiotproof=0)


if 1:
    print(calc.density.nt_sg.shape)
    n = calc.density.nt_sg[0, :, :, 10]
    # 1c4e63
    c0 = np.array([19, 63, 82.0]).reshape((3, 1, 1)) / 255
    c1 = np.array([1.0, 1, 0]).reshape((3, 1, 1))
    a = c0 + n / n.max() * (c1 - c0)
    import pylab as p
    print(a.shape)
    i = p.imshow(a.T, aspect=True)
    i.write_png('ase.png')
    p.axis('off')
    p.savefig('ase2.png', dpi=200)
