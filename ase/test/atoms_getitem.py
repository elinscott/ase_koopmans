from ase.atoms import Atoms


w = Atoms('H2O',
          positions=[[2.264, 0.639, 0.876],
                     [0.792, 0.955, 0.608],
                     [1.347, 0.487, 1.234]],
          cell=[3, 3, 3],
          pbc=True)

try:
    print(w[True, False])
except IndexError:
    pass
else:
    raise

assert(w[0, 1] == w[True, True, False])
assert(w[0, 1] == w[0:2])

