# tests of the dlpoly I/O
from ase_koopmans import io as ase_koopmansIO
from ase_koopmans.io.dlp4 import iread_dlp_history
from io import StringIO

import numpy as np

#Test HISTORY reading
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


def test_dlp():
    cells = []
    cells.append(np.array([[23.01, -0.3943E-01, 0.4612E-01], [-0.9486E-01, 22.98, 0.4551], [0.6568, 0.7694, 19.21]]))
    cells.append(np.array([[22.90, -0.3925E-01, 0.4591E-01], [-0.9443E-01, 22.88, 0.4531], [0.6538, 0.7660, 19.12]]))
    cells.append(np.array([[22.73, -0.3896E-01, 0.4557E-01], [-0.9374E-01, 22.71, 0.4497], [0.6490, 0.7603, 18.98]]))



    traj = ase_koopmansIO.read(fd, format='dlp-history', index=slice(0,None))
    assert len(traj) == 3

    traj = ase_koopmansIO.iread(fd, format='dlp-history', index=slice(0,None))
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

#Test REVCON reading
fd2 = StringIO(u"""                                                                           ch3cl
         2         0         5   103.350212873    
c1               1
    0.1843387826E-03    0.9416060951E-04     1.246412527    
 -0.458762453613E-03  0.115950998393E-02  0.398056206380E-02
 -0.474964352912E-01  0.182081702320      0.229243689875    
cl               2
   -0.2059439421E-03    0.1412072888E-04   -0.5580793078    
 -0.940840256180E-04  0.395951857541E-04 -0.239995080203E-02
 -0.111106651883      0.751717008835E-01  -1.80980665091    
hc               3
   -0.5196456534       -0.9003296761         1.595983262    
 -0.172904226091E-03 -0.532458921776E-03  0.870577192032E-03
 -0.946235115819E-01 -0.277617843971      0.570439906153    
hc               4
     1.040017099        0.7706891154E-04     1.595602116    
  0.842082010780E-03 -0.190324565710E-03  0.148901125710E-02
  0.309680778021     -0.636899601725E-01  0.690354198940    
h1               5
   -0.5195887524        0.9005809880         1.595908521    
 -0.630491688097E-03  0.709696007109E-03  0.401075989715E-02
 -0.564541792647E-01  0.840544009392E-01  0.319768855941    
""")

def test_dlp2():
    mol = ase_koopmansIO.read(fd2, format='dlp4', symbols=['C', 'Cl', 'H', 'H', 'H'])
    assert (mol.get_array('dlp4_labels') == np.array(['1', 'cl', 'hc', 'hc', '1'])).all()
