# tests of the dlpoly I/O
from ase import io as aseIO
from ase.io.dlp4 import iread_dlp_history
from io import StringIO

import numpy as np

fd = StringIO(u"""                                                                             tes
         2         3         2
timestep         1         2         2         3    0.000500
   23.01     -0.3943E-01  0.4612E-01
 -0.9486E-01   22.98      0.4551    
  0.6568      0.7694       19.21    
o1               1   16.000000   -0.730000
  7.9029E+00 -3.7677E+00  4.1862E+00
 -6.6695E-01  1.4565E+00 -1.0286E+01
 -2.9693E+04  4.8692E+04 -3.2588E+05
Ni+               2   16.000000   -0.730000
  7.7028E+00 -3.4396E+00  1.9907E+00
  6.6695E-01 -1.4565E+00  1.0286E+01
  2.9693E+04 -4.8692E+04  3.2588E+05
timestep         2         2         2         3    0.000500
   22.90     -0.3925E-01  0.4591E-01
 -0.9443E-01   22.88      0.4531    
  0.6538      0.7660       19.12    
o1               1   16.000000   -0.730000
  7.9019E+00 -3.7659E+00  4.1735E+00
 -1.5895E+00  2.9698E+00 -2.0415E+01
 -2.9065E+04  4.7608E+04 -3.1855E+05
Ni+               2   16.000000   -0.730000
  7.7038E+00 -3.4414E+00  2.0034E+00
  1.5895E+00 -2.9698E+00  2.0415E+01
  2.9065E+04 -4.7608E+04  3.1855E+05
timestep         3         2         2         3    0.000500
   22.73     -0.3896E-01  0.4557E-01
 -0.9374E-01   22.71      0.4497    
  0.6490      0.7603       18.98    
o1               1   16.000000   -0.730000
  7.9001E+00 -3.7628E+00  4.1528E+00
 -2.4898E+00  4.4453E+00 -3.0289E+01
 -2.8009E+04  4.5827E+04 -3.0655E+05
Ni+               2   16.000000   -0.730000
  7.7056E+00 -3.4445E+00  2.0241E+00
  2.4898E+00 -4.4453E+00  3.0289E+01
  2.8009E+04 -4.5827E+04  3.0655E+05
""")

cells = []
cells.append(np.array([[23.01, -0.3943E-01, 0.4612E-01], [-0.9486E-01, 22.98, 0.4551], [0.6568, 0.7694, 19.21]]))
cells.append(np.array([[22.90, -0.3925E-01, 0.4591E-01], [-0.9443E-01, 22.88, 0.4531], [0.6538, 0.7660, 19.12]]))
cells.append(np.array([[22.73, -0.3896E-01, 0.4557E-01], [-0.9374E-01, 22.71, 0.4497], [0.6490, 0.7603, 18.98]]))



traj = aseIO.read(fd, format='dlp-history', index=slice(0,None))
assert len(traj) == 3

traj = aseIO.iread(fd, format='dlp-history', index=slice(0,None))
for i, frame in enumerate(traj):
    assert len(frame) == 2
    assert frame[0].symbol == 'O'
    assert frame[1].symbol == 'Ni'
    assert np.isclose(frame.get_cell(),cells[i]).all()

symbols = frame.get_chemical_symbols()

traj = iread_dlp_history(fd, symbols)
for i, frame in enumerate(traj):
    assert len(frame) == 2
    assert frame[0].symbol == 'O'
    assert frame[1].symbol == 'Ni'
    assert np.isclose(frame.get_cell(),cells[i]).all()
