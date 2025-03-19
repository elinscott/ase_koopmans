# creates: vasp_si_bandstructure.png
# flake8: noqa
import numpy as np

from ase_koopmans.build import bulk
from ase_koopmans.dft.band_structure import BandStructure

atoms = bulk('Si')

ref = 5.92456665

kpts = np.array([[0.5, 0.25, 0.75],
                 [0.466667, 0.233333, 0.7],
                 [0.433333, 0.216667, 0.65],
                 [0.4, 0.2, 0.6],
                 [0.366667, 0.183333, 0.55],
                 [0.333333, 0.166667, 0.5],
                 [0.3, 0.15, 0.45],
                 [0.266667, 0.133333, 0.4],
                 [0.233333, 0.116667, 0.35],
                 [0.2, 0.1, 0.3],
                 [0.166667, 0.083333, 0.25],
                 [0.133333, 0.066667, 0.2],
                 [0.1, 0.05, 0.15],
                 [0.066667, 0.033333, 0.1],
                 [0.033333, 0.016667, 0.05],
                 [0., 0., 0.],
                 [0.035714, 0., 0.035714],
                 [0.071429, 0., 0.071429],
                 [0.107143, 0., 0.107143],
                 [0.142857, 0., 0.142857],
                 [0.178571, 0., 0.178571],
                 [0.214286, 0., 0.214286],
                 [0.25, 0., 0.25],
                 [0.285714, 0., 0.285714],
                 [0.321429, 0., 0.321429],
                 [0.357143, 0., 0.357143],
                 [0.392857, 0., 0.392857],
                 [0.428571, 0., 0.428571],
                 [0.464286, 0., 0.464286],
                 [0.5, 0., 0.5]])

energies = np.array([[[-1.8719, -1.8719, 1.896, 1.896, 9.9847, 9.9848,
                       10.7179, 10.7179, 16.3164, 16.3164, 18.7643, 18.7644,
                       22.12, 22.12, 24.2139, 24.214, 24.6634, 24.6635,
                       25.647, 25.6471, 29.354, 29.3541, 33.1805, 33.183],
                      [-2.2103, -1.5473, 1.7299, 2.15, 9.7729, 9.8081,
                       10.484, 11.4377, 15.6543, 16.6937, 18.3613, 19.3913,
                       21.486, 22.5026, 23.7336, 23.8698, 25.0892, 25.0928,
                       25.5, 26.0423, 28.835, 29.9382, 31.8389, 33.5773],
                      [-2.6037, -1.1843, 1.6932, 2.443, 9.4486, 9.6039,
                       10.2967, 12.4433, 14.9254, 16.5, 18.5475, 19.9207,
                       20.8098, 22.3632, 22.5643, 24.2574, 25.4175, 25.579,
                       26.0956, 26.217, 28.6439, 30.1375, 30.9884, 33.8185],
                      [-3.0294, -0.7837, 1.7689, 2.7667, 9.1118, 9.4829,
                       10.1514, 13.5061, 14.1964, 15.9525, 19.0523, 20.0908,
                       20.1664, 21.3861, 22.2674, 24.3882, 25.6385, 25.8868,
                       26.3955, 27.0953, 28.9185, 29.1932, 31.5653, 33.7943],
                      [-3.4653, -0.3452, 1.9377, 3.1135, 8.8189, 9.421,
                       10.042, 13.4889, 14.5853, 15.3467, 19.3402, 19.6508,
                       19.8789, 20.443, 22.4422, 23.8594, 25.6162, 26.3879,
                       27.1776, 27.3199, 29.1227, 29.4439, 32.3577, 33.7513],
                      [-3.8933, 0.1309, 2.1829, 3.4731, 8.573, 9.4254,
                       9.9602, 12.8083, 14.7624, 15.6129, 18.7207, 18.9656,
                       19.5516, 20.6621, 22.2717, 23.6605, 25.3157, 26.7113,
                       26.7774, 28.4047, 29.9985, 30.2285, 32.28, 33.2525],
                      [-4.301, 0.6439, 2.4905, 3.833, 8.3732, 9.4968,
                       9.8914, 12.1533, 14.2222, 16.4197, 17.6265, 18.4395,
                       19.3934, 21.227, 21.4938, 24.1161, 24.6403, 26.5225,
                       27.0032, 29.7385, 30.485, 30.8192, 31.6171, 34.122],
                      [-4.6791, 1.1934, 2.8505, 4.1822, 8.2234, 9.5997,
                       9.8702, 11.5362, 13.754, 16.5268, 16.7294, 17.9964,
                       19.6089, 20.7654, 21.6316, 24.0116, 24.7733, 26.3903,
                       27.3109, 29.5881, 30.8624, 31.1446, 32.9522, 34.5147],
                      [-5.0215, 1.7782, 3.2521, 4.5102, 8.1239, 9.581,
                       10.0307, 10.9566, 13.3538, 15.5005, 16.5291, 17.623,
                       19.3788, 21.0737, 21.4625, 23.9596, 25.4904, 26.0844,
                       27.6384, 28.6624, 31.2283, 32.4231, 34.1963, 34.653],
                      [-5.3242, 2.3958, 3.6847, 4.8092, 8.072, 9.4505,
                       10.3284, 10.4101, 12.9888, 14.6112, 16.0752, 17.3176,
                       18.7107, 20.7833, 22.1802, 24.3883, 25.4674, 26.2416,
                       27.9822, 28.2529, 31.4991, 33.2501, 34.6253, 35.0368],
                      [-5.584, 3.0417, 4.1356, 5.0757, 8.0637, 9.2674,
                       9.8997, 10.6249, 12.6042, 13.9822, 15.5114, 17.0853,
                       18.1162, 19.9478, 23.4681, 24.1582, 25.4628, 27.0066,
                       28.2168, 28.3413, 31.7029, 33.2457, 34.1533, 36.0611],
                      [-5.7989, 3.7093, 4.5899, 5.3091, 8.096, 9.0571,
                       9.4335, 10.6789, 12.363, 13.6838, 14.914, 16.9356,
                       17.6552, 19.1182, 23.2204, 24.8241, 26.2215, 27.7893,
                       28.3512, 28.72, 31.84, 32.9072, 33.5178, 36.746],
                      [-5.9674, 4.3838, 5.025, 5.5043, 8.1571, 8.8399,
                       9.015, 10.259, 12.5694, 13.5845, 14.3484, 16.8668,
                       17.3185, 18.3352, 22.3058, 26.1801, 27.0723, 28.5404,
                       28.5538, 29.1176, 31.8815, 32.4588, 32.8549, 36.1559],
                      [-6.0884, 5.0297, 5.4077, 5.6563, 8.2299, 8.6296,
                       8.6612, 9.65, 13.0092, 13.5712, 13.8791, 16.8645,
                       17.0921, 17.6434, 21.5479, 27.4742, 27.8998, 28.7282,
                       29.2679, 29.513, 31.788, 31.996, 32.1979, 35.6287],
                      [-6.1613, 5.5613, 5.6846, 5.7548, 8.2903, 8.4091,
                       8.4278, 9.1627, 13.4162, 13.5657, 13.5826, 16.8969,
                       16.9598, 17.1224, 21.0366, 28.5011, 28.6017, 28.8677,
                       29.856, 29.8939, 31.5295, 31.5761, 31.615, 35.2712],
                      [-6.1857, 5.788, 5.788, 5.788, 8.3138, 8.3138,
                       8.3138, 8.9886, 13.4551, 13.5894, 13.5894, 16.9136,
                       16.9136, 16.9136, 20.8535, 28.8986, 28.8986, 28.8986,
                       30.129, 30.129, 31.3402, 31.3402, 31.3402, 35.1457],
                      [-6.1632, 5.6364, 5.6943, 5.6943, 8.2561, 8.4362,
                       8.4363, 9.1376, 13.3996, 13.5569, 13.6141, 16.956,
                       16.956, 17.0455, 21.0238, 28.5605, 28.5605, 28.947,
                       29.8448, 29.8853, 31.5229, 31.582, 31.5821, 35.3151],
                      [-6.0962, 5.2452, 5.4536, 5.4536, 8.0963, 8.7583,
                       8.7583, 9.5017, 12.9469, 13.6885, 13.8478, 17.0747,
                       17.0747, 17.346, 21.518, 27.7696, 27.7696, 29.0306,
                       29.1391, 29.3657, 31.8584, 32.0218, 32.0218, 35.8194],
                      [-5.9848, 4.7235, 5.1415, 5.1415, 7.867, 9.2049,
                       9.2049, 9.9071, 12.4501, 13.8123, 14.2939, 17.2717,
                       17.2717, 17.6989, 22.2862, 26.8308, 26.8308, 28.2314,
                       28.7248, 29.175, 32.1499, 32.4076, 32.4077, 36.6204],
                      [-5.8292, 4.1392, 4.8081, 4.8081, 7.6025, 9.7251,
                       9.7251, 10.0992, 12.2157, 13.9855, 14.8586, 17.5453,
                       17.5453, 18.0685, 23.256, 25.8576, 25.8576, 27.2453,
                       28.0625, 29.3886, 32.3698, 32.7055, 32.7055, 37.6285],
                      [-5.63, 3.5226, 4.478, 4.478, 7.3293, 9.8663,
                       10.2895, 10.2895, 12.4666, 14.2031, 15.5102, 17.8865,
                       17.8865, 18.4681, 24.3556, 24.8769, 24.8769, 26.2185,
                       27.3801, 29.6856, 32.4636, 32.9117, 32.9118, 37.2868],
                      [-5.3877, 2.8919, 4.1691, 4.1691, 7.0678, 9.3792,
                       10.8851, 10.8851, 13.0349, 14.4747, 16.2214, 18.2921,
                       18.2922, 18.8997, 23.9508, 23.9508, 25.2176, 25.5044,
                       26.7008, 30.096, 32.3961, 33.0677, 33.0678, 36.1881],
                      [-5.103, 2.2538, 3.8875, 3.8875, 6.8315, 8.8348,
                       11.5065, 11.5065, 13.7326, 14.7954, 16.9671, 18.7368,
                       18.7368, 19.3867, 23.1032, 23.1032, 24.2379, 26.0493,
                       26.626, 30.6671, 32.109, 33.0584, 33.0585, 35.1901],
                      [-4.7771, 1.6134, 3.6374, 3.6374, 6.6279, 8.3095,
                       12.1377, 12.1377, 14.475, 15.1559, 17.7119, 19.1372,
                       19.1372, 19.9371, 22.4029, 22.4029, 23.288, 25.4467,
                       27.6348, 31.4424, 31.5389, 32.5669, 32.5669, 34.585],
                      [-4.4111, 0.9752, 3.4224, 3.4224, 6.4645, 7.8328,
                       12.7829, 12.7829, 15.2308, 15.5738, 18.4078, 19.3181,
                       19.3181, 20.5601, 22.0382, 22.0382, 22.3731, 24.9429,
                       28.4734, 30.6975, 31.5772, 31.5772, 32.306, 33.4205],
                      [-4.006, 0.3431, 3.2436, 3.2436, 6.3475, 7.4156,
                       13.4387, 13.4387, 15.9802, 16.0402, 18.9762, 19.0567,
                       19.0567, 21.2593, 21.4901, 22.2082, 22.2082, 24.616,
                       29.1372, 29.6942, 30.3953, 30.3953, 31.679, 33.3345],
                      [-3.5642, -0.2784, 3.1042, 3.1042, 6.2793, 7.064,
                       14.1111, 14.1111, 16.5581, 16.7126, 18.5182, 18.5182,
                       19.3103, 20.6566, 22.0339, 22.8057, 22.8057, 24.5452,
                       28.6214, 29.1794, 29.1794, 29.6229, 30.5584, 32.8451],
                      [-3.0876, -0.8863, 3.0037, 3.0037, 6.2623, 6.7765,
                       14.7728, 14.7728, 17.1201, 17.4077, 17.8504, 17.8504,
                       19.3349, 19.8627, 22.8644, 23.6249, 23.6249, 24.7429,
                       27.5978, 27.9651, 27.9651, 29.0422, 30.4427, 32.3347],
                      [-2.5785, -1.4744, 2.9437, 2.9437, 6.3021, 6.5565,
                       15.41, 15.41, 17.2036, 17.2036, 17.7555, 18.0538,
                       19.0775, 19.1146, 23.7456, 24.5645, 24.5645, 24.9649,
                       26.8166, 26.8167, 26.9285, 27.9975, 30.8961, 31.8679],
                      [-2.0397, -2.0397, 2.9235, 2.9235, 6.399, 6.3991,
                       15.775, 15.775, 16.8298, 16.8298, 18.4119, 18.4119,
                       18.6257, 18.6257, 24.5738, 24.5739, 25.3583, 25.3583,
                       25.9521, 25.9521, 27.1466, 27.1466, 31.3822, 31.3825]]])

# Update to new band structure stuff
lattice = atoms.cell.get_bravais_lattice()
bandpath = lattice.bandpath('WGX', npoints=30)
maxerr = np.abs(bandpath.kpts - kpts).max()
assert maxerr < 1e-5


bs = BandStructure(bandpath,
                   energies=energies,
                   reference=ref)

bs.plot(emin=-13, filename='vasp_si_bandstructure.png')
