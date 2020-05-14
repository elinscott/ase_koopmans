# creates:  diffusion-path.png
from ase.io import read, write
if 1:
    exec(compile(open('diffusion4.py').read(), 'diffusion4.py', 'exec'))
images = [read('mep%d.traj' % i) for i in range(5)]
a = images[0] + images[1] + images[2] + images[3] + images[4]
del a.constraints
a *= (2, 1, 1)
a.set_cell(images[0].get_cell())
write('diffusion-path.pov', a, rotation='-90x',
      transparent=False, run_povray=True)
