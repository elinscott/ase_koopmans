from ase.io.formats import ioformats

traj = ioformats['traj']
print(traj)

outcar = ioformats['vasp-out']
print(outcar)
assert outcar.match_name('OUTCAR')
assert outcar.match_name('something.with.OUTCAR.stuff')


for name, ioformat in ioformats.items():
    print(name)
    print('=' * len(name))
    print(ioformat.full_description())
    print()
