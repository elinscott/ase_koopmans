# creates: C2H4.png
from ase_koopmans.build.molecule import molecule
from ase_koopmans.io import write
from ase_koopmans.io.pov import get_bondpairs, set_high_bondorder_pairs

C2H4 = molecule('C2H4')
r = [{'C': 0.4, 'H': 0.2}[at.symbol] for at in C2H4]
bondpairs = get_bondpairs(C2H4, radius=1.1)
high_bondorder_pairs = {}
# This defines offset, bond order, and bond_offset of the bond between 0 and 1
high_bondorder_pairs[(0, 1)] = ((0, 0, 0), 2, (0.17, 0.17, 0))
bondpairs = set_high_bondorder_pairs(bondpairs, high_bondorder_pairs)

write('C2H4.pov', C2H4, format='pov', run_povray=True,
      canvas_width=200, radii=r, bondatoms=bondpairs, rotation="90y")
