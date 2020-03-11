import numpy as np
from numpy.linalg import norm

from ase.io import read


def test_parse_socketio():
    write_output_socketio()
    traj = read("aims.out", ":", format="aims-output")

    a1, a2 = traj[0], traj[1]
    f1, f2, = a1.get_forces(), a2.get_forces()

    assert a1.positions[0, 0] == 0.01
    assert a2.positions[0, 0] == 0.00

    assert np.allclose(f1[0, 0], -0.119012361951726e00)
    assert np.allclose(f2[0, 0], 0.454895003232345e-03)


def test_run():
    write_output()
    atoms = read("aims.out", format="aims-output")

    # find total energy in aims.out
    key = "| Total energy corrected        :"
    with open("aims.out") as f:
        line = next(l for l in f if key in l)
        ref_energy = float(line.split()[5])

    assert norm(atoms.get_total_energy() - ref_energy) < 1e-12

    # find force in aims.out
    key = "Total atomic forces (unitary forces cleaned) [eV/Ang]:"
    with open("aims.out") as f:
        next(l for l in f if key in l)
        line = next(f)
        ref_force = [float(l) for l in line.split()[2:5]]

    assert norm(atoms.get_forces()[0] - ref_force) < 1e-12

    # find stress in aims.out
    key = "Analytical stress tensor - Symmetrized"
    with open("aims.out") as f:
        next(l for l in f if key in l)
        # scroll to significant lines
        for _ in range(4):
            next(f)
        line = next(f)
        ref_stress = [float(l) for l in line.split()[2:5]]

    assert norm(atoms.get_stress(voigt=False)[0] - ref_stress) < 1e-12

    # find atomic stress in aims.out
    key = "Per atom stress (eV) used for heat flux calculation"
    with open("aims.out") as f:
        next(l for l in f if key in l)
        # scroll to boundary
        next(l for l in f if "-------------" in l)

        line = next(f)
        xx, yy, zz, xy, xz, yz = [float(l) for l in line.split()[2:8]]
        ref_stresses = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]

    assert norm(atoms.get_stresses()[0] - ref_stresses) < 1e-12


def write_output():
    output = "  Basic array size parameters:\n  | Number of species                 :        1\n  | Number of atoms                   :        8\n  | Number of lattice vectors         :        3\n  | Max. basis fn. angular momentum   :        2\n  | Max. atomic/ionic basis occupied n:        3\n  | Max. number of basis fn. types    :        3\n  | Max. radial fns per species/type  :        5\n  | Max. logarithmic grid size        :     1346\n  | Max. radial integration grid size :       42\n  | Max. angular integration grid size:      302\n  | Max. angular grid division number :        8\n  | Radial grid for Hartree potential :     1346\n  | Number of spin channels           :        1\n\n\n  Input geometry:\n  | Unit cell:\n  |        5.42606753        0.00000000        0.00000000\n  |        0.00000000        5.42606753        0.00000000\n  |        0.00000000        0.00000000        5.42606753\n  | Atomic structure:\n  |       Atom                x [A]            y [A]            z [A]\n  |    1: Species Si            0.03431851       -0.09796859        0.09930953\n  |    2: Species Si            5.44231311        2.73920529        2.78205416\n  |    3: Species Si            2.75321969        0.10000784        2.72715717\n  |    4: Species Si            2.73199531        2.68826367       -0.08575931\n  |    5: Species Si            1.34757448        1.42946424        1.25761431\n  |    6: Species Si            1.35486030        4.13154987        4.06589071\n  |    7: Species Si            4.04177845        1.27675199        4.00805480\n  |    8: Species Si            3.99821025        4.01092826        1.42388121\n\n  +-------------------------------------------------------------------+\n  |              Analytical stress tensor - Symmetrized               |\n  |                  Cartesian components [eV/A**3]                   |\n  +-------------------------------------------------------------------+\n  |                x                y                z                |\n  |                                                                   |\n  |  x        -0.01478211      -0.01327277      -0.00355870           |\n  |  y        -0.01327277      -0.01512112      -0.01367280           |\n  |  z        -0.00355870      -0.01367280      -0.01534158           |\n  |                                                                   |\n  |  Pressure:       0.01508160   [eV/A**3]                           |\n  |                                                                   |\n  +-------------------------------------------------------------------+\n\n  ESTIMATED overall HOMO-LUMO gap:      0.21466369 eV between HOMO at k-point 1 and LUMO at k-point 1\n\n  Energy and forces in a compact form:\n  | Total energy uncorrected      :         -0.630943948216411E+05 eV\n  | Total energy corrected        :         -0.630943948568205E+05 eV  <-- do not rely on this value for anything but (periodic) metals\n  | Electronic free energy        :         -0.630943948919999E+05 eV\n  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n  |   1         -0.104637839735875E+01          0.500412824184706E+00         -0.439789552504239E+00\n  |   2         -0.155820611394662E+00         -0.476557335046913E+00         -0.655396151432312E+00\n  |   3         -0.193381405004926E+01         -0.122454085397628E+01         -0.169259060410046E+01\n  |   4          0.404969041951871E-01          0.457139849737633E+00         -0.128445757910440E+00\n  |   5          0.109984435024380E-01         -0.165609149153507E+00          0.114351292468512E+01\n  |   6          0.663029766776301E+00         -0.814079627100908E-01          0.384378715376525E-04\n  |   7          0.213211510059627E+01          0.918575437083381E+00          0.189666102862743E+01\n  |   8          0.289372843732474E+00          0.719871898810707E-01         -0.123990325236629E+00\n\n\n    - Per atom stress (eV) used for heat flux calculation:\n        Atom   | Stress components (1,1), (2,2), (3,3), (1,2), (1,3), (2,3)\n      -------------------------------------------------------------------\n             1 |     0.9843662637E-01   -0.1027274769E+00    0.7237959330E-01   -0.3532042840E+00    0.2563317062E+00   -0.3642257991E+00\n             2 |     0.1244911861E+00   -0.4107147872E-01   -0.1084329966E+00    0.1201650287E+00   -0.1716383020E+00   -0.4669712541E-01\n             3 |    -0.1019986539E+01   -0.7054557814E+00   -0.8410240482E+00   -0.3714228752E+00   -0.4921256188E+00   -0.7970402772E+00\n             4 |    -0.5372048581E+00   -0.2498902919E+00   -0.2260340202E+00   -0.4368600591E+00    0.8622059429E-01    0.9182206824E-01\n             5 |    -0.3268304136E-01   -0.1853638313E+00    0.8046857169E-01   -0.3825550863E+00    0.3088175411E+00   -0.2399437437E+00\n             6 |    -0.2682129292E+00   -0.3832959470E+00   -0.5895171406E+00   -0.8151368635E-02    0.5046578049E-01   -0.6756388823E+00\n             7 |    -0.6970248515E+00   -0.6819450154E+00   -0.9123466446E+00   -0.5254451278E+00   -0.5070403877E+00   -0.6281674944E+00\n             8 |    -0.2933806554E-01   -0.6593089867E-01    0.7360641037E-01   -0.1629233327E+00   -0.9955320981E-01    0.4755870988E+00\n      -------------------------------------------------------------------\n\n\n          Have a nice day.\n------------------------------------------------------------\n"

    with open("aims.out", "w") as f:
        f.write(output)


def write_output_socketio():
    output = """
  MPI-parallelism will be employed.
------------------------------------------------------------
          Invoking FHI-aims ...

          When using FHI-aims, please cite the following reference:

          Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu,
          Ville Havu, Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
          'Ab Initio Molecular Simulations with Numeric Atom-Centered Orbitals',
          Computer Physics Communications 180, 2175-2196 (2009)

          For any questions about FHI-aims, please visit the aimsclub website
          with its forums and wiki. Contributions to both the forums and the
          wiki are warmly encouraged - they are for you, and everyone is welcome there.

------------------------------------------------------------



  Date     :  20191125, Time     :  161252.036
  Time zero on CPU 1             :   0.217210000000000E+00  s.
  Internal wall clock time zero  :           343930372.036  s.

  FHI-aims created a unique identifier for this run for later identification
  aims_uuid : AF6E505A-1838-46FE-A3AE-E2D990E9677F

  Build configuration of the current instance of FHI-aims
  -------------------------------------------------------
  FHI-aims version      : 191023
  Commit number         : 406e2a38d
  CMake host system     : Linux-4.12.14-95.32-default
  CMake version         : 3.5.2
  Fortran compiler      : /mpcdf/soft/SLE_12_SP4/packages/x86_64/intel_parallel_studio/2018.4/compilers_and_libraries_2018.5.274/linux/mpi/bin64/mpiifort (Intel) version 18.0.5.20180823
  Fortran compiler flags: -O3 -ip -fp-model precise
  C compiler            : /mpcdf/soft/SLE_12_SP4/packages/haswell/intel/18.0.5/bin/icc (Intel) version 18.0.5.20180823
  C compiler flags      : -O3
  C++ compiler          : /mpcdf/soft/SLE_12_SP4/packages/haswell/intel/18.0.5/bin/icpc (Intel) version 18.0.5.20180823
  C++ compiler flags    : -O3
  Architecture          : AMD64_AVX
  Using MPI
  Using Scalapack
  Using C files
  Using LibXC
  Using SPGlib
  Using i-PI
  Linking against: /mpcdf/soft/SLE_12_SP4/packages/x86_64/intel_parallel_studio/2018.4/mkl/lib/intel64/libmkl_intel_lp64.so
                   /mpcdf/soft/SLE_12_SP4/packages/x86_64/intel_parallel_studio/2018.4/mkl/lib/intel64/libmkl_sequential.so
                   /mpcdf/soft/SLE_12_SP4/packages/x86_64/intel_parallel_studio/2018.4/mkl/lib/intel64/libmkl_core.so
                   /mpcdf/soft/SLE_12_SP4/packages/x86_64/intel_parallel_studio/2018.4/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so
                   /mpcdf/soft/SLE_12_SP4/packages/x86_64/intel_parallel_studio/2018.4/mkl/lib/intel64/libmkl_scalapack_lp64.so

  Using      192 parallel tasks.
  Task        0 on host dra0334 reporting.
  Task        1 on host dra0334 reporting.
  Task        2 on host dra0334 reporting.
  Task        3 on host dra0334 reporting.
  Task        4 on host dra0334 reporting.
  Task        5 on host dra0334 reporting.
  Task        6 on host dra0334 reporting.
  Task        7 on host dra0334 reporting.
  Task        8 on host dra0334 reporting.
  Task        9 on host dra0334 reporting.
  Task       10 on host dra0334 reporting.
  Task       11 on host dra0334 reporting.
  Task       12 on host dra0334 reporting.
  Task       13 on host dra0334 reporting.
  Task       14 on host dra0334 reporting.
  Task       15 on host dra0334 reporting.
  Task       16 on host dra0334 reporting.
  Task       17 on host dra0334 reporting.
  Task       18 on host dra0334 reporting.
  Task       19 on host dra0334 reporting.
  Task       20 on host dra0334 reporting.
  Task       21 on host dra0334 reporting.
  Task       22 on host dra0334 reporting.
  Task       23 on host dra0334 reporting.
  Task       24 on host dra0334 reporting.
  Task       25 on host dra0334 reporting.
  Task       26 on host dra0334 reporting.
  Task       27 on host dra0334 reporting.
  Task       28 on host dra0334 reporting.
  Task       29 on host dra0334 reporting.
  Task       30 on host dra0334 reporting.
  Task       31 on host dra0334 reporting.
  Task       32 on host dra0337 reporting.
  Task       33 on host dra0337 reporting.
  Task       34 on host dra0337 reporting.
  Task       35 on host dra0337 reporting.
  Task       36 on host dra0337 reporting.
  Task       37 on host dra0337 reporting.
  Task       38 on host dra0337 reporting.
  Task       39 on host dra0337 reporting.
  Task       40 on host dra0337 reporting.
  Task       41 on host dra0337 reporting.
  Task       42 on host dra0337 reporting.
  Task       43 on host dra0337 reporting.
  Task       44 on host dra0337 reporting.
  Task       45 on host dra0337 reporting.
  Task       46 on host dra0337 reporting.
  Task       47 on host dra0337 reporting.
  Task       48 on host dra0337 reporting.
  Task       49 on host dra0337 reporting.
  Task       50 on host dra0337 reporting.
  Task       51 on host dra0337 reporting.
  Task       52 on host dra0337 reporting.
  Task       53 on host dra0337 reporting.
  Task       54 on host dra0337 reporting.
  Task       55 on host dra0337 reporting.
  Task       56 on host dra0337 reporting.
  Task       57 on host dra0337 reporting.
  Task       58 on host dra0337 reporting.
  Task       59 on host dra0337 reporting.
  Task       60 on host dra0337 reporting.
  Task       61 on host dra0337 reporting.
  Task       62 on host dra0337 reporting.
  Task       63 on host dra0337 reporting.
  Task       64 on host dra0339 reporting.
  Task       65 on host dra0339 reporting.
  Task       66 on host dra0339 reporting.
  Task       67 on host dra0339 reporting.
  Task       68 on host dra0339 reporting.
  Task       69 on host dra0339 reporting.
  Task       70 on host dra0339 reporting.
  Task       71 on host dra0339 reporting.
  Task       72 on host dra0339 reporting.
  Task       73 on host dra0339 reporting.
  Task       74 on host dra0339 reporting.
  Task       75 on host dra0339 reporting.
  Task       76 on host dra0339 reporting.
  Task       77 on host dra0339 reporting.
  Task       78 on host dra0339 reporting.
  Task       79 on host dra0339 reporting.
  Task       80 on host dra0339 reporting.
  Task       81 on host dra0339 reporting.
  Task       82 on host dra0339 reporting.
  Task       83 on host dra0339 reporting.
  Task       84 on host dra0339 reporting.
  Task       85 on host dra0339 reporting.
  Task       86 on host dra0339 reporting.
  Task       87 on host dra0339 reporting.
  Task       88 on host dra0339 reporting.
  Task       89 on host dra0339 reporting.
  Task       90 on host dra0339 reporting.
  Task       91 on host dra0339 reporting.
  Task       92 on host dra0339 reporting.
  Task       93 on host dra0339 reporting.
  Task       94 on host dra0339 reporting.
  Task       95 on host dra0339 reporting.
  Task       96 on host dra0340 reporting.
  Task       97 on host dra0340 reporting.
  Task       98 on host dra0340 reporting.
  Task       99 on host dra0340 reporting.
  Task      100 on host dra0340 reporting.
  Task      101 on host dra0340 reporting.
  Task      102 on host dra0340 reporting.
  Task      103 on host dra0340 reporting.
  Task      104 on host dra0340 reporting.
  Task      105 on host dra0340 reporting.
  Task      106 on host dra0340 reporting.
  Task      107 on host dra0340 reporting.
  Task      108 on host dra0340 reporting.
  Task      109 on host dra0340 reporting.
  Task      110 on host dra0340 reporting.
  Task      111 on host dra0340 reporting.
  Task      112 on host dra0340 reporting.
  Task      113 on host dra0340 reporting.
  Task      114 on host dra0340 reporting.
  Task      115 on host dra0340 reporting.
  Task      116 on host dra0340 reporting.
  Task      117 on host dra0340 reporting.
  Task      118 on host dra0340 reporting.
  Task      119 on host dra0340 reporting.
  Task      120 on host dra0340 reporting.
  Task      121 on host dra0340 reporting.
  Task      122 on host dra0340 reporting.
  Task      123 on host dra0340 reporting.
  Task      124 on host dra0340 reporting.
  Task      125 on host dra0340 reporting.
  Task      126 on host dra0340 reporting.
  Task      127 on host dra0340 reporting.
  Task      128 on host dra0341 reporting.
  Task      129 on host dra0341 reporting.
  Task      130 on host dra0341 reporting.
  Task      131 on host dra0341 reporting.
  Task      132 on host dra0341 reporting.
  Task      133 on host dra0341 reporting.
  Task      134 on host dra0341 reporting.
  Task      135 on host dra0341 reporting.
  Task      136 on host dra0341 reporting.
  Task      137 on host dra0341 reporting.
  Task      138 on host dra0341 reporting.
  Task      139 on host dra0341 reporting.
  Task      140 on host dra0341 reporting.
  Task      141 on host dra0341 reporting.
  Task      142 on host dra0341 reporting.
  Task      143 on host dra0341 reporting.
  Task      144 on host dra0341 reporting.
  Task      145 on host dra0341 reporting.
  Task      146 on host dra0341 reporting.
  Task      147 on host dra0341 reporting.
  Task      148 on host dra0341 reporting.
  Task      149 on host dra0341 reporting.
  Task      150 on host dra0341 reporting.
  Task      151 on host dra0341 reporting.
  Task      152 on host dra0341 reporting.
  Task      153 on host dra0341 reporting.
  Task      154 on host dra0341 reporting.
  Task      155 on host dra0341 reporting.
  Task      156 on host dra0341 reporting.
  Task      157 on host dra0341 reporting.
  Task      158 on host dra0341 reporting.
  Task      159 on host dra0341 reporting.
  Task      160 on host dra0345 reporting.
  Task      161 on host dra0345 reporting.
  Task      162 on host dra0345 reporting.
  Task      163 on host dra0345 reporting.
  Task      164 on host dra0345 reporting.
  Task      165 on host dra0345 reporting.
  Task      166 on host dra0345 reporting.
  Task      167 on host dra0345 reporting.
  Task      168 on host dra0345 reporting.
  Task      169 on host dra0345 reporting.
  Task      170 on host dra0345 reporting.
  Task      171 on host dra0345 reporting.
  Task      172 on host dra0345 reporting.
  Task      173 on host dra0345 reporting.
  Task      174 on host dra0345 reporting.
  Task      175 on host dra0345 reporting.
  Task      176 on host dra0345 reporting.
  Task      177 on host dra0345 reporting.
  Task      178 on host dra0345 reporting.
  Task      179 on host dra0345 reporting.
  Task      180 on host dra0345 reporting.
  Task      181 on host dra0345 reporting.
  Task      182 on host dra0345 reporting.
  Task      183 on host dra0345 reporting.
  Task      184 on host dra0345 reporting.
  Task      185 on host dra0345 reporting.
  Task      186 on host dra0345 reporting.
  Task      187 on host dra0345 reporting.
  Task      188 on host dra0345 reporting.
  Task      189 on host dra0345 reporting.
  Task      190 on host dra0345 reporting.
  Task      191 on host dra0345 reporting.

  Performing system and environment tests:
  | Environment variable OMP_NUM_THREADS correctly set to 1.
  | Maximum stacksize for task 0: unlimited
  | Maximum stacksize for task 1: unlimited
  | Maximum stacksize for task 2: unlimited
  | Maximum stacksize for task 3: unlimited
  | Maximum stacksize for task 4: unlimited
  | Maximum stacksize for task 5: unlimited
  | Maximum stacksize for task 6: unlimited
  | Maximum stacksize for task 7: unlimited
  | Maximum stacksize for task 8: unlimited
  | Maximum stacksize for task 9: unlimited
  | Maximum stacksize for task 10: unlimited
  | Maximum stacksize for task 11: unlimited
  | Maximum stacksize for task 12: unlimited
  | Maximum stacksize for task 13: unlimited
  | Maximum stacksize for task 14: unlimited
  | Maximum stacksize for task 15: unlimited
  | Maximum stacksize for task 16: unlimited
  | Maximum stacksize for task 17: unlimited
  | Maximum stacksize for task 18: unlimited
  | Maximum stacksize for task 19: unlimited
  | Maximum stacksize for task 20: unlimited
  | Maximum stacksize for task 21: unlimited
  | Maximum stacksize for task 22: unlimited
  | Maximum stacksize for task 23: unlimited
  | Maximum stacksize for task 24: unlimited
  | Maximum stacksize for task 25: unlimited
  | Maximum stacksize for task 26: unlimited
  | Maximum stacksize for task 27: unlimited
  | Maximum stacksize for task 28: unlimited
  | Maximum stacksize for task 29: unlimited
  | Maximum stacksize for task 30: unlimited
  | Maximum stacksize for task 31: unlimited
  | Maximum stacksize for task 32: unlimited
  | Maximum stacksize for task 33: unlimited
  | Maximum stacksize for task 34: unlimited
  | Maximum stacksize for task 35: unlimited
  | Maximum stacksize for task 36: unlimited
  | Maximum stacksize for task 37: unlimited
  | Maximum stacksize for task 38: unlimited
  | Maximum stacksize for task 39: unlimited
  | Maximum stacksize for task 40: unlimited
  | Maximum stacksize for task 41: unlimited
  | Maximum stacksize for task 42: unlimited
  | Maximum stacksize for task 43: unlimited
  | Maximum stacksize for task 44: unlimited
  | Maximum stacksize for task 45: unlimited
  | Maximum stacksize for task 46: unlimited
  | Maximum stacksize for task 47: unlimited
  | Maximum stacksize for task 48: unlimited
  | Maximum stacksize for task 49: unlimited
  | Maximum stacksize for task 50: unlimited
  | Maximum stacksize for task 51: unlimited
  | Maximum stacksize for task 52: unlimited
  | Maximum stacksize for task 53: unlimited
  | Maximum stacksize for task 54: unlimited
  | Maximum stacksize for task 55: unlimited
  | Maximum stacksize for task 56: unlimited
  | Maximum stacksize for task 57: unlimited
  | Maximum stacksize for task 58: unlimited
  | Maximum stacksize for task 59: unlimited
  | Maximum stacksize for task 60: unlimited
  | Maximum stacksize for task 61: unlimited
  | Maximum stacksize for task 62: unlimited
  | Maximum stacksize for task 63: unlimited
  | Maximum stacksize for task 64: unlimited
  | Maximum stacksize for task 65: unlimited
  | Maximum stacksize for task 66: unlimited
  | Maximum stacksize for task 67: unlimited
  | Maximum stacksize for task 68: unlimited
  | Maximum stacksize for task 69: unlimited
  | Maximum stacksize for task 70: unlimited
  | Maximum stacksize for task 71: unlimited
  | Maximum stacksize for task 72: unlimited
  | Maximum stacksize for task 73: unlimited
  | Maximum stacksize for task 74: unlimited
  | Maximum stacksize for task 75: unlimited
  | Maximum stacksize for task 76: unlimited
  | Maximum stacksize for task 77: unlimited
  | Maximum stacksize for task 78: unlimited
  | Maximum stacksize for task 79: unlimited
  | Maximum stacksize for task 80: unlimited
  | Maximum stacksize for task 81: unlimited
  | Maximum stacksize for task 82: unlimited
  | Maximum stacksize for task 83: unlimited
  | Maximum stacksize for task 84: unlimited
  | Maximum stacksize for task 85: unlimited
  | Maximum stacksize for task 86: unlimited
  | Maximum stacksize for task 87: unlimited
  | Maximum stacksize for task 88: unlimited
  | Maximum stacksize for task 89: unlimited
  | Maximum stacksize for task 90: unlimited
  | Maximum stacksize for task 91: unlimited
  | Maximum stacksize for task 92: unlimited
  | Maximum stacksize for task 93: unlimited
  | Maximum stacksize for task 94: unlimited
  | Maximum stacksize for task 95: unlimited
  | Maximum stacksize for task 96: unlimited
  | Maximum stacksize for task 97: unlimited
  | Maximum stacksize for task 98: unlimited
  | Maximum stacksize for task 99: unlimited
  | Maximum stacksize for task 100: unlimited
  | Maximum stacksize for task 101: unlimited
  | Maximum stacksize for task 102: unlimited
  | Maximum stacksize for task 103: unlimited
  | Maximum stacksize for task 104: unlimited
  | Maximum stacksize for task 105: unlimited
  | Maximum stacksize for task 106: unlimited
  | Maximum stacksize for task 107: unlimited
  | Maximum stacksize for task 108: unlimited
  | Maximum stacksize for task 109: unlimited
  | Maximum stacksize for task 110: unlimited
  | Maximum stacksize for task 111: unlimited
  | Maximum stacksize for task 112: unlimited
  | Maximum stacksize for task 113: unlimited
  | Maximum stacksize for task 114: unlimited
  | Maximum stacksize for task 115: unlimited
  | Maximum stacksize for task 116: unlimited
  | Maximum stacksize for task 117: unlimited
  | Maximum stacksize for task 118: unlimited
  | Maximum stacksize for task 119: unlimited
  | Maximum stacksize for task 120: unlimited
  | Maximum stacksize for task 121: unlimited
  | Maximum stacksize for task 122: unlimited
  | Maximum stacksize for task 123: unlimited
  | Maximum stacksize for task 124: unlimited
  | Maximum stacksize for task 125: unlimited
  | Maximum stacksize for task 126: unlimited
  | Maximum stacksize for task 127: unlimited
  | Maximum stacksize for task 128: unlimited
  | Maximum stacksize for task 129: unlimited
  | Maximum stacksize for task 130: unlimited
  | Maximum stacksize for task 131: unlimited
  | Maximum stacksize for task 132: unlimited
  | Maximum stacksize for task 133: unlimited
  | Maximum stacksize for task 134: unlimited
  | Maximum stacksize for task 135: unlimited
  | Maximum stacksize for task 136: unlimited
  | Maximum stacksize for task 137: unlimited
  | Maximum stacksize for task 138: unlimited
  | Maximum stacksize for task 139: unlimited
  | Maximum stacksize for task 140: unlimited
  | Maximum stacksize for task 141: unlimited
  | Maximum stacksize for task 142: unlimited
  | Maximum stacksize for task 143: unlimited
  | Maximum stacksize for task 144: unlimited
  | Maximum stacksize for task 145: unlimited
  | Maximum stacksize for task 146: unlimited
  | Maximum stacksize for task 147: unlimited
  | Maximum stacksize for task 148: unlimited
  | Maximum stacksize for task 149: unlimited
  | Maximum stacksize for task 150: unlimited
  | Maximum stacksize for task 151: unlimited
  | Maximum stacksize for task 152: unlimited
  | Maximum stacksize for task 153: unlimited
  | Maximum stacksize for task 154: unlimited
  | Maximum stacksize for task 155: unlimited
  | Maximum stacksize for task 156: unlimited
  | Maximum stacksize for task 157: unlimited
  | Maximum stacksize for task 158: unlimited
  | Maximum stacksize for task 159: unlimited
  | Maximum stacksize for task 160: unlimited
  | Maximum stacksize for task 161: unlimited
  | Maximum stacksize for task 162: unlimited
  | Maximum stacksize for task 163: unlimited
  | Maximum stacksize for task 164: unlimited
  | Maximum stacksize for task 165: unlimited
  | Maximum stacksize for task 166: unlimited
  | Maximum stacksize for task 167: unlimited
  | Maximum stacksize for task 168: unlimited
  | Maximum stacksize for task 169: unlimited
  | Maximum stacksize for task 170: unlimited
  | Maximum stacksize for task 171: unlimited
  | Maximum stacksize for task 172: unlimited
  | Maximum stacksize for task 173: unlimited
  | Maximum stacksize for task 174: unlimited
  | Maximum stacksize for task 175: unlimited
  | Maximum stacksize for task 176: unlimited
  | Maximum stacksize for task 177: unlimited
  | Maximum stacksize for task 178: unlimited
  | Maximum stacksize for task 179: unlimited
  | Maximum stacksize for task 180: unlimited
  | Maximum stacksize for task 181: unlimited
  | Maximum stacksize for task 182: unlimited
  | Maximum stacksize for task 183: unlimited
  | Maximum stacksize for task 184: unlimited
  | Maximum stacksize for task 185: unlimited
  | Maximum stacksize for task 186: unlimited
  | Maximum stacksize for task 187: unlimited
  | Maximum stacksize for task 188: unlimited
  | Maximum stacksize for task 189: unlimited
  | Maximum stacksize for task 190: unlimited
  | Maximum stacksize for task 191: unlimited
  | Current stacksize for task 0: unlimited
  | Current stacksize for task 1: unlimited
  | Current stacksize for task 2: unlimited
  | Current stacksize for task 3: unlimited
  | Current stacksize for task 4: unlimited
  | Current stacksize for task 5: unlimited
  | Current stacksize for task 6: unlimited
  | Current stacksize for task 7: unlimited
  | Current stacksize for task 8: unlimited
  | Current stacksize for task 9: unlimited
  | Current stacksize for task 10: unlimited
  | Current stacksize for task 11: unlimited
  | Current stacksize for task 12: unlimited
  | Current stacksize for task 13: unlimited
  | Current stacksize for task 14: unlimited
  | Current stacksize for task 15: unlimited
  | Current stacksize for task 16: unlimited
  | Current stacksize for task 17: unlimited
  | Current stacksize for task 18: unlimited
  | Current stacksize for task 19: unlimited
  | Current stacksize for task 20: unlimited
  | Current stacksize for task 21: unlimited
  | Current stacksize for task 22: unlimited
  | Current stacksize for task 23: unlimited
  | Current stacksize for task 24: unlimited
  | Current stacksize for task 25: unlimited
  | Current stacksize for task 26: unlimited
  | Current stacksize for task 27: unlimited
  | Current stacksize for task 28: unlimited
  | Current stacksize for task 29: unlimited
  | Current stacksize for task 30: unlimited
  | Current stacksize for task 31: unlimited
  | Current stacksize for task 32: unlimited
  | Current stacksize for task 33: unlimited
  | Current stacksize for task 34: unlimited
  | Current stacksize for task 35: unlimited
  | Current stacksize for task 36: unlimited
  | Current stacksize for task 37: unlimited
  | Current stacksize for task 38: unlimited
  | Current stacksize for task 39: unlimited
  | Current stacksize for task 40: unlimited
  | Current stacksize for task 41: unlimited
  | Current stacksize for task 42: unlimited
  | Current stacksize for task 43: unlimited
  | Current stacksize for task 44: unlimited
  | Current stacksize for task 45: unlimited
  | Current stacksize for task 46: unlimited
  | Current stacksize for task 47: unlimited
  | Current stacksize for task 48: unlimited
  | Current stacksize for task 49: unlimited
  | Current stacksize for task 50: unlimited
  | Current stacksize for task 51: unlimited
  | Current stacksize for task 52: unlimited
  | Current stacksize for task 53: unlimited
  | Current stacksize for task 54: unlimited
  | Current stacksize for task 55: unlimited
  | Current stacksize for task 56: unlimited
  | Current stacksize for task 57: unlimited
  | Current stacksize for task 58: unlimited
  | Current stacksize for task 59: unlimited
  | Current stacksize for task 60: unlimited
  | Current stacksize for task 61: unlimited
  | Current stacksize for task 62: unlimited
  | Current stacksize for task 63: unlimited
  | Current stacksize for task 64: unlimited
  | Current stacksize for task 65: unlimited
  | Current stacksize for task 66: unlimited
  | Current stacksize for task 67: unlimited
  | Current stacksize for task 68: unlimited
  | Current stacksize for task 69: unlimited
  | Current stacksize for task 70: unlimited
  | Current stacksize for task 71: unlimited
  | Current stacksize for task 72: unlimited
  | Current stacksize for task 73: unlimited
  | Current stacksize for task 74: unlimited
  | Current stacksize for task 75: unlimited
  | Current stacksize for task 76: unlimited
  | Current stacksize for task 77: unlimited
  | Current stacksize for task 78: unlimited
  | Current stacksize for task 79: unlimited
  | Current stacksize for task 80: unlimited
  | Current stacksize for task 81: unlimited
  | Current stacksize for task 82: unlimited
  | Current stacksize for task 83: unlimited
  | Current stacksize for task 84: unlimited
  | Current stacksize for task 85: unlimited
  | Current stacksize for task 86: unlimited
  | Current stacksize for task 87: unlimited
  | Current stacksize for task 88: unlimited
  | Current stacksize for task 89: unlimited
  | Current stacksize for task 90: unlimited
  | Current stacksize for task 91: unlimited
  | Current stacksize for task 92: unlimited
  | Current stacksize for task 93: unlimited
  | Current stacksize for task 94: unlimited
  | Current stacksize for task 95: unlimited
  | Current stacksize for task 96: unlimited
  | Current stacksize for task 97: unlimited
  | Current stacksize for task 98: unlimited
  | Current stacksize for task 99: unlimited
  | Current stacksize for task 100: unlimited
  | Current stacksize for task 101: unlimited
  | Current stacksize for task 102: unlimited
  | Current stacksize for task 103: unlimited
  | Current stacksize for task 104: unlimited
  | Current stacksize for task 105: unlimited
  | Current stacksize for task 106: unlimited
  | Current stacksize for task 107: unlimited
  | Current stacksize for task 108: unlimited
  | Current stacksize for task 109: unlimited
  | Current stacksize for task 110: unlimited
  | Current stacksize for task 111: unlimited
  | Current stacksize for task 112: unlimited
  | Current stacksize for task 113: unlimited
  | Current stacksize for task 114: unlimited
  | Current stacksize for task 115: unlimited
  | Current stacksize for task 116: unlimited
  | Current stacksize for task 117: unlimited
  | Current stacksize for task 118: unlimited
  | Current stacksize for task 119: unlimited
  | Current stacksize for task 120: unlimited
  | Current stacksize for task 121: unlimited
  | Current stacksize for task 122: unlimited
  | Current stacksize for task 123: unlimited
  | Current stacksize for task 124: unlimited
  | Current stacksize for task 125: unlimited
  | Current stacksize for task 126: unlimited
  | Current stacksize for task 127: unlimited
  | Current stacksize for task 128: unlimited
  | Current stacksize for task 129: unlimited
  | Current stacksize for task 130: unlimited
  | Current stacksize for task 131: unlimited
  | Current stacksize for task 132: unlimited
  | Current stacksize for task 133: unlimited
  | Current stacksize for task 134: unlimited
  | Current stacksize for task 135: unlimited
  | Current stacksize for task 136: unlimited
  | Current stacksize for task 137: unlimited
  | Current stacksize for task 138: unlimited
  | Current stacksize for task 139: unlimited
  | Current stacksize for task 140: unlimited
  | Current stacksize for task 141: unlimited
  | Current stacksize for task 142: unlimited
  | Current stacksize for task 143: unlimited
  | Current stacksize for task 144: unlimited
  | Current stacksize for task 145: unlimited
  | Current stacksize for task 146: unlimited
  | Current stacksize for task 147: unlimited
  | Current stacksize for task 148: unlimited
  | Current stacksize for task 149: unlimited
  | Current stacksize for task 150: unlimited
  | Current stacksize for task 151: unlimited
  | Current stacksize for task 152: unlimited
  | Current stacksize for task 153: unlimited
  | Current stacksize for task 154: unlimited
  | Current stacksize for task 155: unlimited
  | Current stacksize for task 156: unlimited
  | Current stacksize for task 157: unlimited
  | Current stacksize for task 158: unlimited
  | Current stacksize for task 159: unlimited
  | Current stacksize for task 160: unlimited
  | Current stacksize for task 161: unlimited
  | Current stacksize for task 162: unlimited
  | Current stacksize for task 163: unlimited
  | Current stacksize for task 164: unlimited
  | Current stacksize for task 165: unlimited
  | Current stacksize for task 166: unlimited
  | Current stacksize for task 167: unlimited
  | Current stacksize for task 168: unlimited
  | Current stacksize for task 169: unlimited
  | Current stacksize for task 170: unlimited
  | Current stacksize for task 171: unlimited
  | Current stacksize for task 172: unlimited
  | Current stacksize for task 173: unlimited
  | Current stacksize for task 174: unlimited
  | Current stacksize for task 175: unlimited
  | Current stacksize for task 176: unlimited
  | Current stacksize for task 177: unlimited
  | Current stacksize for task 178: unlimited
  | Current stacksize for task 179: unlimited
  | Current stacksize for task 180: unlimited
  | Current stacksize for task 181: unlimited
  | Current stacksize for task 182: unlimited
  | Current stacksize for task 183: unlimited
  | Current stacksize for task 184: unlimited
  | Current stacksize for task 185: unlimited
  | Current stacksize for task 186: unlimited
  | Current stacksize for task 187: unlimited
  | Current stacksize for task 188: unlimited
  | Current stacksize for task 189: unlimited
  | Current stacksize for task 190: unlimited
  | Current stacksize for task 191: unlimited
  | Checking for scalapack...
  | Testing pdtran()...
  | All pdtran() tests passed.

  Obtaining array dimensions for all initial allocations:

  -----------------------------------------------------------------------
  Parsing control.in (first pass over file, find array dimensions only).
  The contents of control.in will be repeated verbatim below
  unless switched off by setting 'verbatim_writeout .false.' .
  in the first line of control.in .
  -----------------------------------------------------------------------

  #===============================================================================
  # FHI-aims file: ./control.in
  # Created using the Atomic Simulation Environment (ASE)
  # Mon Nov 25 16:12:48 2019
  #
  # List of parameters used to initialize the calculator:
  #     xc : pbesol
  #     charge_mix_param : 0.3
  #     sc_accuracy_rho : 1e-06
  #     relativistic : atomic_zora scalar
  #     output_level : MD_light
  #     compute_forces : True
  #     k_grid : [2, 2, 2]
  #     use_pimd_wrapper : ('localhost', 12345)
  #     species_dir : /u/fknoop/FHI-aims/aimsfiles/species_defaults/light
  #===============================================================================
  xc                                 pbesol
  charge_mix_param                   0.3
  sc_accuracy_rho                    1e-06
  relativistic                       atomic_zora scalar
  output_level                       MD_light
  compute_forces                     .true.
  k_grid                             2 2 2
  use_pimd_wrapper                   localhost 12345
  #===============================================================================

  ################################################################################
  #
  #  FHI-aims code project
  #  VB, Fritz-Haber Institut, 2009
  #
  #  Suggested "light" defaults for Mg atom (to be pasted into control.in file)
  #  Be sure to double-check any results obtained with these settings for post-processing,
  #  e.g., with the "tight" defaults and larger basis sets.
  #
  ################################################################################
    species        Mg
  #     global species definitions
      nucleus             12
      mass                24.3050
  #
      l_hartree           4
  #
      cut_pot             4.0          1.5  1.0
      basis_dep_cutoff    1e-4
  #
      radial_base         40 5.5
      radial_multiplier   1
      angular_grids       specified
        division   0.7029   50
        division   0.9689  110
        division   1.1879  194
        division   1.3129  302
  #      division   1.4867  434
  #      division   1.6018  590
  #      division   1.8611  770
  #      division   1.9576  974
  #      division   2.2261 1202
  #      outer_grid   974
        outer_grid   302
  ################################################################################
  #
  #  Definition of "minimal" basis
  #
  ################################################################################
  #     valence basis states
      valence      3  s   2.
      valence      2  p   6.
  #     ion occupancy
      ion_occ      2  s   2.
      ion_occ      2  p   6.
  ################################################################################
  #
  #  Suggested additional basis functions. For production calculations,
  #  uncomment them one after another (the most important basis functions are
  #  listed first).
  #
  #  Constructed for dimers: 2.125 A, 2.375 A, 2.875 A, 3.375 A, 4.5 A
  #
  ################################################################################
  #  "First tier" - improvements: -230.76 meV to -21.94 meV
       hydro 2 p 1.5
       ionic 3 d auto
       hydro 3 s 2.4
  #  "Second tier" - improvements: -5.43 meV to -1.64 meV
  #     hydro 4 f 4.3
  #     hydro 2 p 3.4
  #     hydro 4 s 11.2
  #     hydro 3 d 6.2
  #  "Third tier" - improvements: -0.92 meV to -0.22 meV
  #     hydro 2 s 0.6
  #     hydro 3 p 4.8
  #     hydro 4 f 7.4
  #     hydro 5 g 6.6
  #     hydro 2 p 1.6
  #     hydro 3 d 1.8
  #  "Fourth tier" - improvements: -0.09 meV to -0.05 meV
  #     hydro 4 p 0.45
  #     hydro 5 g 10.4
  #     hydro 2 s 12.4
  #     hydro 4 d 1.7
  ################################################################################
  #
  #  FHI-aims code project
  #  VB, Fritz-Haber Institut, 2009
  #
  #  Suggested "light" defaults for O atom (to be pasted into control.in file)
  #  Be sure to double-check any results obtained with these settings for post-processing,
  #  e.g., with the "tight" defaults and larger basis sets.
  #
  ################################################################################
    species        O
  #     global species definitions
      nucleus             8
      mass                15.9994
  #
      l_hartree           4
  #
      cut_pot             3.5  1.5  1.0
      basis_dep_cutoff    1e-4
  #
      radial_base         36 5.0
      radial_multiplier   1
       angular_grids specified
        division   0.2659   50
        division   0.4451  110
        division   0.6052  194
        division   0.7543  302
  #      division   0.8014  434
  #      division   0.8507  590
  #      division   0.8762  770
  #      division   0.9023  974
  #      division   1.2339 1202
  #      outer_grid 974
        outer_grid 302
  ################################################################################
  #
  #  Definition of "minimal" basis
  #
  ################################################################################
  #     valence basis states
      valence      2  s   2.
      valence      2  p   4.
  #     ion occupancy
      ion_occ      2  s   1.
      ion_occ      2  p   3.
  ################################################################################
  #
  #  Suggested additional basis functions. For production calculations,
  #  uncomment them one after another (the most important basis functions are
  #  listed first).
  #
  #  Constructed for dimers: 1.0 A, 1.208 A, 1.5 A, 2.0 A, 3.0 A
  #
  ################################################################################
  #  "First tier" - improvements: -699.05 meV to -159.38 meV
       hydro 2 p 1.8
       hydro 3 d 7.6
       hydro 3 s 6.4
  #  "Second tier" - improvements: -49.91 meV to -5.39 meV
  #     hydro 4 f 11.6
  #     hydro 3 p 6.2
  #     hydro 3 d 5.6
  #     hydro 5 g 17.6
  #     hydro 1 s 0.75
  #  "Third tier" - improvements: -2.83 meV to -0.50 meV
  #     ionic 2 p auto
  #     hydro 4 f 10.8
  #     hydro 4 d 4.7
  #     hydro 2 s 6.8
  #  "Fourth tier" - improvements: -0.40 meV to -0.12 meV
  #     hydro 3 p 5
  #     hydro 3 s 3.3
  #     hydro 5 g 15.6
  #     hydro 4 f 17.6
  #     hydro 4 d 14
  # Further basis functions - -0.08 meV and below
  #     hydro 3 s 2.1
  #     hydro 4 d 11.6
  #     hydro 3 p 16
  #     hydro 2 s 17.2

  -----------------------------------------------------------------------
  Completed first pass over input file control.in .
  -----------------------------------------------------------------------


  -----------------------------------------------------------------------
  Parsing geometry.in (first pass over file, find array dimensions only).
  The contents of geometry.in will be repeated verbatim below
  unless switched off by setting 'verbatim_writeout .false.' .
  in the first line of geometry.in .
  -----------------------------------------------------------------------

  #=======================================================
  # FHI-aims file: ./geometry.in
  # Created using the Atomic Simulation Environment (ASE)
  # Mon Nov 25 16:12:48 2019
  #=======================================================
  lattice_vector 12.6788144400000000 0.0000000000000002 0.0000000000000000
  lattice_vector 0.0000000000000002 12.6788144400000000 0.0000000000000000
  lattice_vector -0.0000000000000002 -0.0000000000000002 12.6788144400000000
  atom 0.0100000000000000 0.0000000000000000 0.0000000000000000 Mg
  atom 12.6788144400000000 2.1131357400000002 2.1131357399999997 Mg
  atom 12.6788144400000000 4.2262714799999994 4.2262714799999994 Mg
  atom 12.6788144400000000 6.3394072200000000 6.3394072200000000 Mg
  atom 12.6788144400000000 8.4525429599999988 8.4525429599999988 Mg
  atom 12.6788144400000000 10.5656787000000012 10.5656787000000012 Mg
  atom 2.1131357399999997 -0.0000000000000000 2.1131357399999997 Mg
  atom 2.1131357399999997 2.1131357399999997 4.2262714799999994 Mg
  atom 2.1131357399999997 4.2262714799999994 6.3394072200000000 Mg
  atom 2.1131357399999997 6.3394072200000000 8.4525429599999988 Mg
  atom 2.1131357399999997 8.4525429599999988 10.5656786999999994 Mg
  atom 2.1131357399999988 10.5656787000000012 0.0000000000000000 Mg
  atom 4.2262714799999994 -0.0000000000000000 4.2262714799999994 Mg
  atom 4.2262714799999994 2.1131357399999997 6.3394072200000000 Mg
  atom 4.2262714799999994 4.2262714799999994 8.4525429599999988 Mg
  atom 4.2262714799999994 6.3394072200000000 10.5656786999999994 Mg
  atom 4.2262714799999994 8.4525429599999988 0.0000000000000000 Mg
  atom 4.2262714799999994 10.5656787000000012 2.1131357400000010 Mg
  atom 6.3394072200000000 0.0000000000000000 6.3394072200000000 Mg
  atom 6.3394072200000000 2.1131357399999997 8.4525429599999988 Mg
  atom 6.3394072200000000 4.2262714799999994 10.5656786999999994 Mg
  atom 6.3394072200000000 6.3394072200000000 0.0000000000000000 Mg
  atom 6.3394072199999991 8.4525429599999988 2.1131357399999979 Mg
  atom 6.3394072199999991 10.5656787000000012 4.2262714800000021 Mg
  atom 8.4525429599999988 -0.0000000000000000 8.4525429599999988 Mg
  atom 8.4525429599999988 2.1131357399999997 10.5656786999999994 Mg
  atom 8.4525429599999988 4.2262714799999994 0.0000000000000000 Mg
  atom 8.4525429599999988 6.3394072200000000 2.1131357399999979 Mg
  atom 8.4525429599999988 8.4525429599999988 4.2262714799999994 Mg
  atom 8.4525429599999988 10.5656787000000012 6.3394072200000000 Mg
  atom 10.5656787000000012 0.0000000000000000 10.5656787000000012 Mg
  atom 10.5656787000000012 2.1131357399999993 0.0000000000000000 Mg
  atom 10.5656787000000012 4.2262714799999994 2.1131357400000010 Mg
  atom 10.5656787000000012 6.3394072200000000 4.2262714800000021 Mg
  atom 10.5656787000000012 8.4525429599999988 6.3394072200000000 Mg
  atom 10.5656787000000012 10.5656787000000012 8.4525429600000006 Mg
  atom 2.1131357399999997 2.1131357399999997 0.0000000000000000 Mg
  atom 2.1131357399999997 4.2262714799999994 2.1131357399999997 Mg
  atom 2.1131357399999997 6.3394072200000000 4.2262714799999994 Mg
  atom 2.1131357399999997 8.4525429599999988 6.3394072200000000 Mg
  atom 2.1131357399999993 10.5656786999999994 8.4525429599999988 Mg
  atom 2.1131357399999997 -0.0000000000000002 10.5656787000000012 Mg
  atom 4.2262714799999994 2.1131357399999997 2.1131357399999997 Mg
  atom 4.2262714799999994 4.2262714799999994 4.2262714799999994 Mg
  atom 4.2262714799999994 6.3394072200000000 6.3394072200000000 Mg
  atom 4.2262714799999994 8.4525429599999988 8.4525429599999988 Mg
  atom 4.2262714799999994 10.5656786999999994 10.5656786999999994 Mg
  atom 4.2262714799999994 0.0000000000000001 0.0000000000000000 Mg
  atom 6.3394072200000000 2.1131357399999997 4.2262714799999994 Mg
  atom 6.3394072200000000 4.2262714799999994 6.3394072200000000 Mg
  atom 6.3394072200000000 6.3394072200000000 8.4525429599999988 Mg
  atom 6.3394072199999982 8.4525429599999988 10.5656786999999994 Mg
  atom 6.3394072199999982 10.5656786999999994 0.0000000000000000 Mg
  atom 6.3394072199999982 0.0000000000000001 2.1131357400000010 Mg
  atom 8.4525429599999988 2.1131357399999997 6.3394072200000000 Mg
  atom 8.4525429599999988 4.2262714799999994 8.4525429599999988 Mg
  atom 8.4525429599999988 6.3394072200000000 10.5656786999999994 Mg
  atom 8.4525429599999988 8.4525429599999988 0.0000000000000000 Mg
  atom 8.4525429599999988 10.5656786999999994 2.1131357399999979 Mg
  atom 8.4525429599999988 0.0000000000000001 4.2262714800000021 Mg
  atom 10.5656786999999994 2.1131357399999993 8.4525429599999988 Mg
  atom 10.5656786999999994 4.2262714799999994 10.5656786999999994 Mg
  atom 10.5656786999999994 6.3394072200000000 0.0000000000000000 Mg
  atom 10.5656786999999994 8.4525429599999988 2.1131357399999979 Mg
  atom 10.5656786999999994 10.5656786999999994 4.2262714799999994 Mg
  atom 10.5656786999999994 0.0000000000000001 6.3394072200000000 Mg
  atom -0.0000000000000002 2.1131357399999997 10.5656787000000012 Mg
  atom 0.0000000000000001 4.2262714799999994 0.0000000000000000 Mg
  atom 0.0000000000000001 6.3394072200000000 2.1131357400000010 Mg
  atom 0.0000000000000001 8.4525429599999988 4.2262714800000021 Mg
  atom 0.0000000000000001 10.5656786999999994 6.3394072200000000 Mg
  atom -0.0000000000000002 -0.0000000000000002 8.4525429600000006 Mg
  atom 4.2262714799999994 4.2262714799999994 0.0000000000000000 Mg
  atom 4.2262714799999994 6.3394072200000000 2.1131357399999997 Mg
  atom 4.2262714799999994 8.4525429599999988 4.2262714799999994 Mg
  atom 4.2262714799999994 10.5656786999999994 6.3394072200000000 Mg
  atom 4.2262714799999994 -0.0000000000000001 8.4525429599999988 Mg
  atom 4.2262714799999994 2.1131357400000015 10.5656787000000012 Mg
  atom 6.3394072200000000 4.2262714799999994 2.1131357399999997 Mg
  atom 6.3394072200000000 6.3394072200000000 4.2262714799999994 Mg
  atom 6.3394072199999991 8.4525429599999988 6.3394072200000000 Mg
  atom 6.3394072199999991 10.5656786999999994 8.4525429599999988 Mg
  atom 6.3394072199999991 -0.0000000000000001 10.5656786999999994 Mg
  atom 6.3394072199999982 2.1131357400000010 0.0000000000000000 Mg
  atom 8.4525429599999988 4.2262714799999994 4.2262714799999994 Mg
  atom 8.4525429599999988 6.3394072200000000 6.3394072200000000 Mg
  atom 8.4525429599999988 8.4525429599999988 8.4525429599999988 Mg
  atom 8.4525429599999988 10.5656786999999994 10.5656786999999994 Mg
  atom 8.4525429599999988 0.0000000000000002 0.0000000000000000 Mg
  atom 8.4525429599999988 2.1131357400000010 2.1131357400000010 Mg
  atom 10.5656786999999994 4.2262714799999994 6.3394072200000000 Mg
  atom 10.5656786999999994 6.3394072200000000 8.4525429599999988 Mg
  atom 10.5656786999999994 8.4525429599999988 10.5656786999999994 Mg
  atom 10.5656786999999994 10.5656786999999994 0.0000000000000000 Mg
  atom 10.5656786999999994 0.0000000000000002 2.1131357399999979 Mg
  atom 10.5656786999999994 2.1131357400000006 4.2262714800000021 Mg
  atom -0.0000000000000001 4.2262714799999994 8.4525429599999988 Mg
  atom -0.0000000000000001 6.3394072200000000 10.5656786999999994 Mg
  atom 0.0000000000000002 8.4525429599999988 0.0000000000000000 Mg
  atom 0.0000000000000002 10.5656786999999994 2.1131357399999979 Mg
  atom -0.0000000000000001 -0.0000000000000001 4.2262714799999994 Mg
  atom -0.0000000000000001 2.1131357400000010 6.3394072200000000 Mg
  atom 2.1131357400000015 4.2262714799999994 10.5656787000000012 Mg
  atom 2.1131357400000010 6.3394072200000000 0.0000000000000000 Mg
  atom 2.1131357400000010 8.4525429599999988 2.1131357400000010 Mg
  atom 2.1131357400000006 10.5656786999999994 4.2262714800000021 Mg
  atom 2.1131357400000010 -0.0000000000000001 6.3394072200000000 Mg
  atom 2.1131357400000010 2.1131357400000010 8.4525429600000006 Mg
  atom 2.1131357399999997 2.1131357399999997 2.1131357399999997 O
  atom 2.1131357399999997 4.2262714799999994 4.2262714799999994 O
  atom 2.1131357399999997 6.3394072200000000 6.3394072200000000 O
  atom 2.1131357399999997 8.4525429599999988 8.4525429599999988 O
  atom 2.1131357399999997 10.5656787000000012 10.5656787000000012 O
  atom 2.1131357399999993 12.6788144399999982 12.6788144399999982 O
  atom 4.2262714799999994 2.1131357399999997 4.2262714799999994 O
  atom 4.2262714799999994 4.2262714799999994 6.3394072200000000 O
  atom 4.2262714799999994 6.3394072200000000 8.4525429600000006 O
  atom 4.2262714799999994 8.4525429599999988 10.5656786999999994 O
  atom 4.2262714799999994 10.5656787000000012 0.0000000000000000 O
  atom 4.2262714799999994 12.6788144399999982 2.1131357399999979 O
  atom 6.3394072200000000 2.1131357399999997 6.3394072200000000 O
  atom 6.3394072200000000 4.2262714799999994 8.4525429600000006 O
  atom 6.3394072200000000 6.3394072200000000 10.5656787000000012 O
  atom 6.3394072199999991 8.4525429599999988 0.0000000000000000 O
  atom 6.3394072199999991 10.5656787000000012 2.1131357400000010 O
  atom 6.3394072199999991 12.6788144399999982 4.2262714799999994 O
  atom 8.4525429599999988 2.1131357399999997 8.4525429599999988 O
  atom 8.4525429599999988 4.2262714799999994 10.5656786999999994 O
  atom 8.4525429599999988 6.3394072200000000 0.0000000000000000 O
  atom 8.4525429599999988 8.4525429599999988 2.1131357399999979 O
  atom 8.4525429599999988 10.5656787000000012 4.2262714799999994 O
  atom 8.4525429599999988 12.6788144399999982 6.3394072199999973 O
  atom 10.5656787000000012 2.1131357399999997 10.5656787000000012 O
  atom 10.5656787000000012 4.2262714799999994 0.0000000000000000 O
  atom 10.5656787000000012 6.3394072200000000 2.1131357400000010 O
  atom 10.5656787000000012 8.4525429599999988 4.2262714799999994 O
  atom 10.5656787000000012 10.5656787000000012 6.3394072200000000 O
  atom 10.5656787000000012 12.6788144399999982 8.4525429599999988 O
  atom 12.6788144399999982 2.1131357399999997 12.6788144399999982 O
  atom 12.6788144399999982 4.2262714799999994 2.1131357399999979 O
  atom 12.6788144399999982 6.3394072200000000 4.2262714799999994 O
  atom 12.6788144399999982 8.4525429599999988 6.3394072199999973 O
  atom 12.6788144399999982 10.5656787000000012 8.4525429599999988 O
  atom 12.6788144399999982 12.6788144399999982 10.5656786999999959 O
  atom 4.2262714799999994 4.2262714799999994 2.1131357399999997 O
  atom 4.2262714799999994 6.3394072200000000 4.2262714799999994 O
  atom 4.2262714799999994 8.4525429599999988 6.3394072200000000 O
  atom 4.2262714799999994 10.5656786999999994 8.4525429599999988 O
  atom 4.2262714799999994 -0.0000000000000001 10.5656787000000012 O
  atom 4.2262714799999994 2.1131357399999975 12.6788144399999982 O
  atom 6.3394072200000000 4.2262714799999994 4.2262714799999994 O
  atom 6.3394072200000000 6.3394072200000000 6.3394072200000000 O
  atom 6.3394072199999991 8.4525429599999988 8.4525429600000006 O
  atom 6.3394072199999991 10.5656786999999994 10.5656786999999994 O
  atom 6.3394072199999991 0.0000000000000001 0.0000000000000000 O
  atom 6.3394072199999991 2.1131357399999979 2.1131357399999979 O
  atom 8.4525429599999988 4.2262714799999994 6.3394072200000000 O
  atom 8.4525429599999988 6.3394072200000000 8.4525429600000006 O
  atom 8.4525429599999988 8.4525429599999988 10.5656787000000012 O
  atom 8.4525429599999988 10.5656786999999994 0.0000000000000000 O
  atom 8.4525429599999988 0.0000000000000001 2.1131357400000010 O
  atom 8.4525429599999988 2.1131357399999979 4.2262714799999994 O
  atom 10.5656786999999994 4.2262714799999994 8.4525429599999988 O
  atom 10.5656786999999994 6.3394072200000000 10.5656786999999994 O
  atom 10.5656786999999994 8.4525429599999988 0.0000000000000000 O
  atom 10.5656786999999994 10.5656786999999994 2.1131357399999979 O
  atom 10.5656786999999994 0.0000000000000001 4.2262714799999994 O
  atom 10.5656786999999994 2.1131357399999975 6.3394072199999973 O
  atom -0.0000000000000001 4.2262714799999994 10.5656787000000012 O
  atom 0.0000000000000001 6.3394072200000000 0.0000000000000000 O
  atom 0.0000000000000001 8.4525429599999988 2.1131357400000010 O
  atom 0.0000000000000001 10.5656786999999994 4.2262714799999994 O
  atom -0.0000000000000001 -0.0000000000000001 6.3394072200000000 O
  atom -0.0000000000000001 2.1131357399999979 8.4525429599999988 O
  atom 2.1131357399999975 4.2262714799999994 12.6788144399999982 O
  atom 2.1131357399999979 6.3394072200000000 2.1131357399999979 O
  atom 2.1131357399999979 8.4525429599999988 4.2262714799999994 O
  atom 2.1131357399999975 10.5656786999999994 6.3394072199999973 O
  atom 2.1131357399999979 -0.0000000000000001 8.4525429599999988 O
  atom 2.1131357399999979 2.1131357399999979 10.5656786999999959 O
  atom 6.3394072200000000 6.3394072200000000 2.1131357399999997 O
  atom 6.3394072200000000 8.4525429600000006 4.2262714799999994 O
  atom 6.3394072200000000 10.5656787000000012 6.3394072200000000 O
  atom 6.3394072200000000 -0.0000000000000000 8.4525429599999988 O
  atom 6.3394072200000000 2.1131357400000015 10.5656787000000012 O
  atom 6.3394072199999991 4.2262714799999994 12.6788144399999982 O
  atom 8.4525429600000006 6.3394072200000000 4.2262714799999994 O
  atom 8.4525429600000006 8.4525429600000006 6.3394072200000000 O
  atom 8.4525429599999988 10.5656787000000012 8.4525429600000006 O
  atom 8.4525429599999988 -0.0000000000000000 10.5656786999999994 O
  atom 8.4525429599999988 2.1131357400000010 0.0000000000000000 O
  atom 8.4525429599999988 4.2262714799999994 2.1131357399999979 O
  atom 10.5656787000000012 6.3394072200000000 6.3394072200000000 O
  atom 10.5656787000000012 8.4525429600000006 8.4525429600000006 O
  atom 10.5656787000000012 10.5656787000000012 10.5656787000000012 O
  atom 10.5656786999999994 0.0000000000000002 0.0000000000000000 O
  atom 10.5656786999999994 2.1131357400000006 2.1131357400000010 O
  atom 10.5656786999999994 4.2262714799999994 4.2262714799999994 O
  atom -0.0000000000000000 6.3394072200000000 8.4525429599999988 O
  atom -0.0000000000000000 8.4525429600000006 10.5656786999999994 O
  atom 0.0000000000000002 10.5656787000000012 0.0000000000000000 O
  atom -0.0000000000000000 -0.0000000000000000 2.1131357399999979 O
  atom -0.0000000000000000 2.1131357400000010 4.2262714799999994 O
  atom -0.0000000000000000 4.2262714799999994 6.3394072199999973 O
  atom 2.1131357400000015 6.3394072200000000 10.5656787000000012 O
  atom 2.1131357400000010 8.4525429600000006 0.0000000000000000 O
  atom 2.1131357400000006 10.5656787000000012 2.1131357400000010 O
  atom 2.1131357400000010 -0.0000000000000000 4.2262714799999994 O
  atom 2.1131357400000010 2.1131357400000010 6.3394072200000000 O
  atom 2.1131357400000010 4.2262714799999994 8.4525429599999988 O
  atom 4.2262714799999994 6.3394072200000000 12.6788144399999982 O
  atom 4.2262714799999994 8.4525429600000006 2.1131357399999979 O
  atom 4.2262714799999994 10.5656787000000012 4.2262714799999994 O
  atom 4.2262714799999994 -0.0000000000000000 6.3394072199999973 O
  atom 4.2262714799999994 2.1131357400000010 8.4525429599999988 O
  atom 4.2262714799999994 4.2262714799999994 10.5656786999999959 O

  -----------------------------------------------------------------------
  Completed first pass over input file geometry.in .
  -----------------------------------------------------------------------


  Basic array size parameters:
  | Number of species                 :        2
  | Number of atoms                   :      216
  | Number of lattice vectors         :        3
  | Max. basis fn. angular momentum   :        2
  | Max. atomic/ionic basis occupied n:        3
  | Max. number of basis fn. types    :        3
  | Max. radial fns per species/type  :        4
  | Max. logarithmic grid size        :     1334
  | Max. radial integration grid size :       40
  | Max. angular integration grid size:      302
  | Max. angular grid division number :        8
  | Radial grid for Hartree potential :     1334
  | Number of spin channels           :        1

------------------------------------------------------------
          Reading file control.in.
------------------------------------------------------------
  XC: Using PBEsol gradient-corrected functionals.
  Charge density mixing - mixing parameter:     0.3000
  Convergence accuracy of self-consistent charge density:  0.1000E-05
  Scalar relativistic treatment of kinetic energy: on-site free-atom approximation to ZORA.
  Requested output level: MD_light
  Forces will be computed.
  Found k-point grid:         2         2         2
  Using external wrapper (i-PI) for performing (path integral) molecular dynamics
  **Attention: initial geometry.in file will be ignored!

  Reading configuration options for species Mg                  .
  | Found nuclear charge :  12.0000
  | Found atomic mass :    24.3050000000000      amu
  | Found l_max for Hartree potential  :   4
  | Found cutoff potl. onset [A], width [A], scale factor :    4.00000    1.50000    1.00000
  | Threshold for basis-dependent cutoff potential is   0.100000E-03
  | Found data for basic radial integration grid :    40 points, outermost radius =    5.500 A
  | Found multiplier for basic radial grid :   1
  | Found angular grid specification: user-specified.
  | Specified grid contains     5 separate shells.
  | Check grid settings after all constraints further below.
  | Found free-atom valence shell :  3 s   2.000
  | Found free-atom valence shell :  2 p   6.000
  | Found free-ion valence shell :  2 s   2.000
  | Found free-ion valence shell :  2 p   6.000
  | Found hydrogenic basis function :  2 p   1.500
  | Found ionic basis function :  3 d , default cutoff radius.
  | Found hydrogenic basis function :  3 s   2.400
  Species Mg                  : Missing cutoff potential type.
  Defaulting to exp(1/x)/(1-x)^2 type cutoff potential.
  Species Mg: No 'logarithmic' tag. Using default grid for free atom:
  | Default logarithmic grid data [bohr] : 0.1000E-03 0.1000E+03 0.1012E+01
  | Will include ionic basis functions of  2.0-fold positive Mg                   ion.
  Species Mg: On-site basis accuracy parameter (for Gram-Schmidt orthonormalisation) not specified.
  Using default value basis_acc =  0.1000000E-03.
  Species Mg                  : Using default innermost maximum threshold i_radial=  2 for radial functions.
  Species Mg                  : Default cutoff onset for free atom density etc. : 0.40000000E+01 AA.
  Species Mg                  : Basic radial grid will be enhanced according to radial_multiplier =   1, to contain    40 grid points.

  Reading configuration options for species O                   .
  | Found nuclear charge :   8.0000
  | Found atomic mass :    15.9994000000000      amu
  | Found l_max for Hartree potential  :   4
  | Found cutoff potl. onset [A], width [A], scale factor :    3.50000    1.50000    1.00000
  | Threshold for basis-dependent cutoff potential is   0.100000E-03
  | Found data for basic radial integration grid :    36 points, outermost radius =    5.000 A
  | Found multiplier for basic radial grid :   1
  | Found angular grid specification: user-specified.
  | Specified grid contains     5 separate shells.
  | Check grid settings after all constraints further below.
  | Found free-atom valence shell :  2 s   2.000
  | Found free-atom valence shell :  2 p   4.000
  | Found free-ion valence shell :  2 s   1.000
  | Found free-ion valence shell :  2 p   3.000
  | Found hydrogenic basis function :  2 p   1.800
  | Found hydrogenic basis function :  3 d   7.600
  | Found hydrogenic basis function :  3 s   6.400
  Species O                   : Missing cutoff potential type.
  Defaulting to exp(1/x)/(1-x)^2 type cutoff potential.
  Species O : No 'logarithmic' tag. Using default grid for free atom:
  | Default logarithmic grid data [bohr] : 0.1000E-03 0.1000E+03 0.1012E+01
  Species O : On-site basis accuracy parameter (for Gram-Schmidt orthonormalisation) not specified.
  Using default value basis_acc =  0.1000000E-03.
  Species O                   : Using default innermost maximum threshold i_radial=  2 for radial functions.
  Species O                   : Default cutoff onset for free atom density etc. : 0.35000000E+01 AA.
  Species O                   : Basic radial grid will be enhanced according to radial_multiplier =   1, to contain    36 grid points.

  Finished reading input file 'control.in'. Consistency checks are next.

  MPI_IN_PLACE appears to work with this MPI implementation.
  | Keeping use_mpi_in_place .true. (see manual).
  Target number of points in a grid batch is not set. Defaulting to  100
  Method for grid partitioning is not set. Defaulting to parallel hash+maxmin partitioning.
  Batch size limit is not set. Defaulting to    200
  By default, will store active basis functions for each batch.
  If in need of memory, prune_basis_once .false. can be used to disable this option.
  communication_type for Hartree potential was not specified.
  Defaulting to calc_hartree .
  Defaulting to Pulay charge density mixer.
  Pulay mixer: Number of relevant iterations not set.
  Defaulting to    8 iterations.
  Pulay mixer: Number of initial linear mixing iterations not set.
  Defaulting to    0 iterations.
  Work space size for distributed Hartree potential not set.
  Defaulting to   0.200000E+03 MB.
  Algorithm-dependent basis array size parameters:
  | n_max_pulay                         :        8
  Maximum number of self-consistency iterations not provided.
  Presetting  1000 iterations.
  Presetting      1001 iterations before the initial mixing cycle
  is restarted anyway using the sc_init_iter criterion / keyword.
  Presetting a factor      1.000 between actual scf density residual
  and density convergence criterion sc_accuracy_rho below which sc_init_iter
  takes no effect.
  Geometry relaxation not requested: no relaxation will be performed.
  Force calculation: scf convergence accuracy of forces not set.
  Defaulting to 'sc_accuracy_forces not_checked'.
  Handling of forces: Unphysical translation and rotation will be removed from forces.
  No accuracy limit for integral partition fn. given. Defaulting to  0.1000E-14.
  No threshold value for u(r) in integrations given. Defaulting to  0.1000E-05.
  No occupation type (smearing scheme) given. Defaulting to Gaussian broadening, width =  0.1000E-01 eV.
  The width will be adjusted in iteration number     2 of the first full s.c.f. cycle only.
  S.C.F. convergence parameters will be adjusted in iteration number     2 of the first full s.c.f. cycle only.
  No accuracy for occupation numbers given. Defaulting to  0.1000E-12.
  No threshold value for occupation numbers given. Defaulting to  0.0000E+00.
  No accuracy for fermi level given. Defaulting to  0.1000E-19.
  Maximum # of iterations to find E_F not set. Defaulting to  200.
  Preferred method for the eigenvalue solver ('KS_method') not specified in 'control.in'.
  Calling BLACS routine to test compilation state
  Since ScaLAPACK support is enabled, defaulting to ELPA (via ELSI).
  Will not use alltoall communication since running on < 1024 CPUs.
  Threshold for basis singularities not set.
  Default threshold for basis singularities:  0.1000E-04
  partition_type (choice of integration weights) for integrals was not specified.
  | Using a version of the partition function of Stratmann and coworkers ('stratmann_sparse').
  | At each grid point, the set of atoms used to build the partition table is smoothly restricted to
  | only those atoms whose free-atom density would be non-zero at that grid point.
  Partitioning for Hartree potential was not defined. Using partition_type for integrals.
  | Adjusted default value of keyword multip_moments_threshold to:       0.10000000E-11
  | This value may affect high angular momentum components of the Hartree potential in periodic systems.
  Spin handling was not defined in control.in. Defaulting to unpolarized case.
  Angular momentum expansion for Kerker preconditioner not set explicitly.
  | Using default value of   0
  No explicit requirement for turning off preconditioner.
  | By default, it will be turned off when the charge convergence reaches
  | sc_accuracy_rho  =   0.100000E-05
  No special mixing parameter while Kerker preconditioner is on.
  Using default: charge_mix_param =     0.3000.
  No q(lm)/r^(l+1) cutoff set for long-range Hartree potential.
  | Using default value of  0.100000E-09 .
  | Verify using the multipole_threshold keyword.
  Defaulting to new monopole extrapolation.
  Density update method: automatic selection selected.
  Using density matrix based charge density update.
  Using density matrix based charge density update.
  Using packed matrix style: index .
  Defaulting to use time-reversal symmetry for k-point grid.
  Charge integration errors on the 3D integration grid will be compensated
  by explicit normalization and distribution of residual charges.
  Use the "compensate_multipole_errors" flag to change this behaviour.
------------------------------------------------------------


------------------------------------------------------------
          Reading geometry description geometry.in.
------------------------------------------------------------
  Input structure read successfully.
  The structure contains      216 atoms,  and a total of       2160.000 electrons.

  Input geometry:
  | Unit cell:
  |       12.67881444        0.00000000        0.00000000
  |        0.00000000       12.67881444        0.00000000
  |       -0.00000000       -0.00000000       12.67881444
  | Atomic structure:
  |       Atom                x [A]            y [A]            z [A]
  |    1: Species Mg            0.01000000        0.00000000        0.00000000
  |    2: Species Mg           12.67881444        2.11313574        2.11313574
  |    3: Species Mg           12.67881444        4.22627148        4.22627148
  |    4: Species Mg           12.67881444        6.33940722        6.33940722
  |    5: Species Mg           12.67881444        8.45254296        8.45254296
  |    6: Species Mg           12.67881444       10.56567870       10.56567870
  |    7: Species Mg            2.11313574        0.00000000        2.11313574
  |    8: Species Mg            2.11313574        2.11313574        4.22627148
  |    9: Species Mg            2.11313574        4.22627148        6.33940722
  |   10: Species Mg            2.11313574        6.33940722        8.45254296
  |   11: Species Mg            2.11313574        8.45254296       10.56567870
  |   12: Species Mg            2.11313574       10.56567870        0.00000000
  |   13: Species Mg            4.22627148        0.00000000        4.22627148
  |   14: Species Mg            4.22627148        2.11313574        6.33940722
  |   15: Species Mg            4.22627148        4.22627148        8.45254296
  |   16: Species Mg            4.22627148        6.33940722       10.56567870
  |   17: Species Mg            4.22627148        8.45254296        0.00000000
  |   18: Species Mg            4.22627148       10.56567870        2.11313574
  |   19: Species Mg            6.33940722        0.00000000        6.33940722
  |   20: Species Mg            6.33940722        2.11313574        8.45254296
  |   21: Species Mg            6.33940722        4.22627148       10.56567870
  |   22: Species Mg            6.33940722        6.33940722        0.00000000
  |   23: Species Mg            6.33940722        8.45254296        2.11313574
  |   24: Species Mg            6.33940722       10.56567870        4.22627148
  |   25: Species Mg            8.45254296        0.00000000        8.45254296
  |   26: Species Mg            8.45254296        2.11313574       10.56567870
  |   27: Species Mg            8.45254296        4.22627148        0.00000000
  |   28: Species Mg            8.45254296        6.33940722        2.11313574
  |   29: Species Mg            8.45254296        8.45254296        4.22627148
  |   30: Species Mg            8.45254296       10.56567870        6.33940722
  |   31: Species Mg           10.56567870        0.00000000       10.56567870
  |   32: Species Mg           10.56567870        2.11313574        0.00000000
  |   33: Species Mg           10.56567870        4.22627148        2.11313574
  |   34: Species Mg           10.56567870        6.33940722        4.22627148
  |   35: Species Mg           10.56567870        8.45254296        6.33940722
  |   36: Species Mg           10.56567870       10.56567870        8.45254296
  |   37: Species Mg            2.11313574        2.11313574        0.00000000
  |   38: Species Mg            2.11313574        4.22627148        2.11313574
  |   39: Species Mg            2.11313574        6.33940722        4.22627148
  |   40: Species Mg            2.11313574        8.45254296        6.33940722
  |   41: Species Mg            2.11313574       10.56567870        8.45254296
  |   42: Species Mg            2.11313574       -0.00000000       10.56567870
  |   43: Species Mg            4.22627148        2.11313574        2.11313574
  |   44: Species Mg            4.22627148        4.22627148        4.22627148
  |   45: Species Mg            4.22627148        6.33940722        6.33940722
  |   46: Species Mg            4.22627148        8.45254296        8.45254296
  |   47: Species Mg            4.22627148       10.56567870       10.56567870
  |   48: Species Mg            4.22627148        0.00000000        0.00000000
  |   49: Species Mg            6.33940722        2.11313574        4.22627148
  |   50: Species Mg            6.33940722        4.22627148        6.33940722
  |   51: Species Mg            6.33940722        6.33940722        8.45254296
  |   52: Species Mg            6.33940722        8.45254296       10.56567870
  |   53: Species Mg            6.33940722       10.56567870        0.00000000
  |   54: Species Mg            6.33940722        0.00000000        2.11313574
  |   55: Species Mg            8.45254296        2.11313574        6.33940722
  |   56: Species Mg            8.45254296        4.22627148        8.45254296
  |   57: Species Mg            8.45254296        6.33940722       10.56567870
  |   58: Species Mg            8.45254296        8.45254296        0.00000000
  |   59: Species Mg            8.45254296       10.56567870        2.11313574
  |   60: Species Mg            8.45254296        0.00000000        4.22627148
  |   61: Species Mg           10.56567870        2.11313574        8.45254296
  |   62: Species Mg           10.56567870        4.22627148       10.56567870
  |   63: Species Mg           10.56567870        6.33940722        0.00000000
  |   64: Species Mg           10.56567870        8.45254296        2.11313574
  |   65: Species Mg           10.56567870       10.56567870        4.22627148
  |   66: Species Mg           10.56567870        0.00000000        6.33940722
  |   67: Species Mg           -0.00000000        2.11313574       10.56567870
  |   68: Species Mg            0.00000000        4.22627148        0.00000000
  |   69: Species Mg            0.00000000        6.33940722        2.11313574
  |   70: Species Mg            0.00000000        8.45254296        4.22627148
  |   71: Species Mg            0.00000000       10.56567870        6.33940722
  |   72: Species Mg           -0.00000000       -0.00000000        8.45254296
  |   73: Species Mg            4.22627148        4.22627148        0.00000000
  |   74: Species Mg            4.22627148        6.33940722        2.11313574
  |   75: Species Mg            4.22627148        8.45254296        4.22627148
  |   76: Species Mg            4.22627148       10.56567870        6.33940722
  |   77: Species Mg            4.22627148       -0.00000000        8.45254296
  |   78: Species Mg            4.22627148        2.11313574       10.56567870
  |   79: Species Mg            6.33940722        4.22627148        2.11313574
  |   80: Species Mg            6.33940722        6.33940722        4.22627148
  |   81: Species Mg            6.33940722        8.45254296        6.33940722
  |   82: Species Mg            6.33940722       10.56567870        8.45254296
  |   83: Species Mg            6.33940722       -0.00000000       10.56567870
  |   84: Species Mg            6.33940722        2.11313574        0.00000000
  |   85: Species Mg            8.45254296        4.22627148        4.22627148
  |   86: Species Mg            8.45254296        6.33940722        6.33940722
  |   87: Species Mg            8.45254296        8.45254296        8.45254296
  |   88: Species Mg            8.45254296       10.56567870       10.56567870
  |   89: Species Mg            8.45254296        0.00000000        0.00000000
  |   90: Species Mg            8.45254296        2.11313574        2.11313574
  |   91: Species Mg           10.56567870        4.22627148        6.33940722
  |   92: Species Mg           10.56567870        6.33940722        8.45254296
  |   93: Species Mg           10.56567870        8.45254296       10.56567870
  |   94: Species Mg           10.56567870       10.56567870        0.00000000
  |   95: Species Mg           10.56567870        0.00000000        2.11313574
  |   96: Species Mg           10.56567870        2.11313574        4.22627148
  |   97: Species Mg           -0.00000000        4.22627148        8.45254296
  |   98: Species Mg           -0.00000000        6.33940722       10.56567870
  |   99: Species Mg            0.00000000        8.45254296        0.00000000
  |  100: Species Mg            0.00000000       10.56567870        2.11313574
  |  101: Species Mg           -0.00000000       -0.00000000        4.22627148
  |  102: Species Mg           -0.00000000        2.11313574        6.33940722
  |  103: Species Mg            2.11313574        4.22627148       10.56567870
  |  104: Species Mg            2.11313574        6.33940722        0.00000000
  |  105: Species Mg            2.11313574        8.45254296        2.11313574
  |  106: Species Mg            2.11313574       10.56567870        4.22627148
  |  107: Species Mg            2.11313574       -0.00000000        6.33940722
  |  108: Species Mg            2.11313574        2.11313574        8.45254296
  |  109: Species O             2.11313574        2.11313574        2.11313574
  |  110: Species O             2.11313574        4.22627148        4.22627148
  |  111: Species O             2.11313574        6.33940722        6.33940722
  |  112: Species O             2.11313574        8.45254296        8.45254296
  |  113: Species O             2.11313574       10.56567870       10.56567870
  |  114: Species O             2.11313574       12.67881444       12.67881444
  |  115: Species O             4.22627148        2.11313574        4.22627148
  |  116: Species O             4.22627148        4.22627148        6.33940722
  |  117: Species O             4.22627148        6.33940722        8.45254296
  |  118: Species O             4.22627148        8.45254296       10.56567870
  |  119: Species O             4.22627148       10.56567870        0.00000000
  |  120: Species O             4.22627148       12.67881444        2.11313574
  |  121: Species O             6.33940722        2.11313574        6.33940722
  |  122: Species O             6.33940722        4.22627148        8.45254296
  |  123: Species O             6.33940722        6.33940722       10.56567870
  |  124: Species O             6.33940722        8.45254296        0.00000000
  |  125: Species O             6.33940722       10.56567870        2.11313574
  |  126: Species O             6.33940722       12.67881444        4.22627148
  |  127: Species O             8.45254296        2.11313574        8.45254296
  |  128: Species O             8.45254296        4.22627148       10.56567870
  |  129: Species O             8.45254296        6.33940722        0.00000000
  |  130: Species O             8.45254296        8.45254296        2.11313574
  |  131: Species O             8.45254296       10.56567870        4.22627148
  |  132: Species O             8.45254296       12.67881444        6.33940722
  |  133: Species O            10.56567870        2.11313574       10.56567870
  |  134: Species O            10.56567870        4.22627148        0.00000000
  |  135: Species O            10.56567870        6.33940722        2.11313574
  |  136: Species O            10.56567870        8.45254296        4.22627148
  |  137: Species O            10.56567870       10.56567870        6.33940722
  |  138: Species O            10.56567870       12.67881444        8.45254296
  |  139: Species O            12.67881444        2.11313574       12.67881444
  |  140: Species O            12.67881444        4.22627148        2.11313574
  |  141: Species O            12.67881444        6.33940722        4.22627148
  |  142: Species O            12.67881444        8.45254296        6.33940722
  |  143: Species O            12.67881444       10.56567870        8.45254296
  |  144: Species O            12.67881444       12.67881444       10.56567870
  |  145: Species O             4.22627148        4.22627148        2.11313574
  |  146: Species O             4.22627148        6.33940722        4.22627148
  |  147: Species O             4.22627148        8.45254296        6.33940722
  |  148: Species O             4.22627148       10.56567870        8.45254296
  |  149: Species O             4.22627148       -0.00000000       10.56567870
  |  150: Species O             4.22627148        2.11313574       12.67881444
  |  151: Species O             6.33940722        4.22627148        4.22627148
  |  152: Species O             6.33940722        6.33940722        6.33940722
  |  153: Species O             6.33940722        8.45254296        8.45254296
  |  154: Species O             6.33940722       10.56567870       10.56567870
  |  155: Species O             6.33940722        0.00000000        0.00000000
  |  156: Species O             6.33940722        2.11313574        2.11313574
  |  157: Species O             8.45254296        4.22627148        6.33940722
  |  158: Species O             8.45254296        6.33940722        8.45254296
  |  159: Species O             8.45254296        8.45254296       10.56567870
  |  160: Species O             8.45254296       10.56567870        0.00000000
  |  161: Species O             8.45254296        0.00000000        2.11313574
  |  162: Species O             8.45254296        2.11313574        4.22627148
  |  163: Species O            10.56567870        4.22627148        8.45254296
  |  164: Species O            10.56567870        6.33940722       10.56567870
  |  165: Species O            10.56567870        8.45254296        0.00000000
  |  166: Species O            10.56567870       10.56567870        2.11313574
  |  167: Species O            10.56567870        0.00000000        4.22627148
  |  168: Species O            10.56567870        2.11313574        6.33940722
  |  169: Species O            -0.00000000        4.22627148       10.56567870
  |  170: Species O             0.00000000        6.33940722        0.00000000
  |  171: Species O             0.00000000        8.45254296        2.11313574
  |  172: Species O             0.00000000       10.56567870        4.22627148
  |  173: Species O            -0.00000000       -0.00000000        6.33940722
  |  174: Species O            -0.00000000        2.11313574        8.45254296
  |  175: Species O             2.11313574        4.22627148       12.67881444
  |  176: Species O             2.11313574        6.33940722        2.11313574
  |  177: Species O             2.11313574        8.45254296        4.22627148
  |  178: Species O             2.11313574       10.56567870        6.33940722
  |  179: Species O             2.11313574       -0.00000000        8.45254296
  |  180: Species O             2.11313574        2.11313574       10.56567870
  |  181: Species O             6.33940722        6.33940722        2.11313574
  |  182: Species O             6.33940722        8.45254296        4.22627148
  |  183: Species O             6.33940722       10.56567870        6.33940722
  |  184: Species O             6.33940722        0.00000000        8.45254296
  |  185: Species O             6.33940722        2.11313574       10.56567870
  |  186: Species O             6.33940722        4.22627148       12.67881444
  |  187: Species O             8.45254296        6.33940722        4.22627148
  |  188: Species O             8.45254296        8.45254296        6.33940722
  |  189: Species O             8.45254296       10.56567870        8.45254296
  |  190: Species O             8.45254296        0.00000000       10.56567870
  |  191: Species O             8.45254296        2.11313574        0.00000000
  |  192: Species O             8.45254296        4.22627148        2.11313574
  |  193: Species O            10.56567870        6.33940722        6.33940722
  |  194: Species O            10.56567870        8.45254296        8.45254296
  |  195: Species O            10.56567870       10.56567870       10.56567870
  |  196: Species O            10.56567870        0.00000000        0.00000000
  |  197: Species O            10.56567870        2.11313574        2.11313574
  |  198: Species O            10.56567870        4.22627148        4.22627148
  |  199: Species O             0.00000000        6.33940722        8.45254296
  |  200: Species O             0.00000000        8.45254296       10.56567870
  |  201: Species O             0.00000000       10.56567870        0.00000000
  |  202: Species O             0.00000000        0.00000000        2.11313574
  |  203: Species O             0.00000000        2.11313574        4.22627148
  |  204: Species O             0.00000000        4.22627148        6.33940722
  |  205: Species O             2.11313574        6.33940722       10.56567870
  |  206: Species O             2.11313574        8.45254296        0.00000000
  |  207: Species O             2.11313574       10.56567870        2.11313574
  |  208: Species O             2.11313574        0.00000000        4.22627148
  |  209: Species O             2.11313574        2.11313574        6.33940722
  |  210: Species O             2.11313574        4.22627148        8.45254296
  |  211: Species O             4.22627148        6.33940722       12.67881444
  |  212: Species O             4.22627148        8.45254296        2.11313574
  |  213: Species O             4.22627148       10.56567870        4.22627148
  |  214: Species O             4.22627148        0.00000000        6.33940722
  |  215: Species O             4.22627148        2.11313574        8.45254296
  |  216: Species O             4.22627148        4.22627148       10.56567870

  Lattice parameters for 3D lattice (in Angstroms) :    12.678814   12.678814   12.678814
  Angle(s) between unit vectors (in degrees)       :    90.000000   90.000000   90.000000

  | The smallest distance between any two atoms is         2.10313574 AA.
  | The first atom of this pair is atom number                      1 .
  | The second atom of this pair is atom number                   114 .
  | Wigner-Seitz cell of the first atom image           0     0     0 .
  | (The Wigner-Seitz cell of the second atom is 0 0 0  by definition.)

  Symmetry information by spglib:
  | Precision set to  0.1E-04
  | Number of Operations  : 8
  | Space group           : 99
  | International         : P4mm
  | Schoenflies           : C4v^1

  Quantities derived from the lattice vectors:
  | Reciprocal lattice vector 1:  0.495566  0.000000  0.000000
  | Reciprocal lattice vector 2:  0.000000  0.495566  0.000000
  | Reciprocal lattice vector 3:  0.000000  0.000000  0.495566
  | Unit cell volume                               :   0.203815E+04  A^3

  Range separation radius for Ewald summation (hartree_convergence_parameter):      2.90574234 bohr.

  Number of empty states per atom not set in control.in - providing a guess from actual geometry.
  | Total number of empty states used during s.c.f. cycle:      648
  If you use a very high smearing, use empty_states (per atom!) in control.in to increase this value.

  Structure-dependent array size parameters:
  | Maximum number of distinct radial functions  :       13
  | Maximum number of basis functions            :     3132
  | Number of Kohn-Sham states (occupied + empty):     1728
------------------------------------------------------------

************************** W A R N I N G *******************************
* You are using the PIMD wrapper. Specifications and positions         *
* in geometry.in will be IGNORED - all is received from the wrapper.   *
* Please make sure species are declared in the same order in           *
* geometry.in and wrapper input.                                       *
************************************************************************


------------------------------------------------------------
          Preparing all fixed parts of the calculation.
------------------------------------------------------------
  Determining machine precision:
    2.225073858507201E-308
  Setting up grids for atomic and cluster calculations.

  Creating wave function, potential, and density for free atoms.

  Species: Mg

  List of occupied orbitals and eigenvalues:
    n    l              occ      energy [Ha]    energy [eV]
    1    0           2.0000       -46.359170     -1261.4972
    2    0           2.0000        -2.927498       -79.6613
    3    0           2.0000        -0.167982        -4.5710
    2    1           6.0000        -1.706030       -46.4234


  Species: O

  List of occupied orbitals and eigenvalues:
    n    l              occ      energy [Ha]    energy [eV]
    1    0           2.0000       -18.926989      -515.0296
    2    0           2.0000        -0.880247       -23.9527
    2    1           4.0000        -0.331514        -9.0210

  Creating fixed part of basis set: Ionic, confined, hydrogenic.

  Mg                   ion:

  List of free ionic orbitals and eigenvalues:
    n    l      energy [Ha]    energy [eV]
    1    0       -47.122987     -1282.2817
    2    0        -3.674099       -99.9773
    2    1        -2.451441       -66.7071


  List of ionic basis orbitals and eigenvalues:
    n    l      energy [Ha]    energy [eV]    outer radius [A]
    3    2        -0.253797        -6.9062       5.100062


  Mg                   hydrogenic:

  List of hydrogenic basis orbitals:
    n    l      effective z      eigenvalue [eV]  inner max. [A]     outer max. [A]     outer radius [A]
    2    1         1.500000        -7.5955           1.395712           1.395712           5.100062
    3    0         2.400000        -8.0111           0.162319           2.734035           5.100062


  O                    hydrogenic:

  List of hydrogenic basis orbitals:
    n    l      effective z      eigenvalue [eV]  inner max. [A]     outer max. [A]     outer radius [A]
    2    1         1.800000       -10.9749           1.164242           1.164242           4.578029
    3    2         7.600000       -87.3180           0.624125           0.624125           3.251020
    3    0         6.400000       -61.9207           0.061167           1.081902           4.001998


  Adding cutoff potential to free-atom effective potential.
  Creating atomic-like basis functions for current effective potential.

  Species Mg                  :

  List of atomic basis orbitals and eigenvalues:
    n    l      energy [Ha]    energy [eV]    outer radius [A]
    1    0       -46.359170     -1261.4972       0.921046
    2    0        -2.927498       -79.6613       3.621735
    3    0        -0.167982        -4.5710       5.100062
    2    1        -1.706030       -46.4234       4.404175


  Species O                   :

  List of atomic basis orbitals and eigenvalues:
    n    l      energy [Ha]    energy [eV]    outer radius [A]
    1    0       -18.926989      -515.0296       1.415765
    2    0        -0.880247       -23.9527       4.413171
    2    1        -0.331514        -9.0210       4.522403

  Assembling full basis from fixed parts.
  | Species Mg :   atomic orbital   1 s accepted.
  | Species Mg :   atomic orbital   2 s accepted.
  | Species Mg :   atomic orbital   3 s accepted.
  | Species Mg :    hydro orbital   3 s accepted.
  | Species Mg :   atomic orbital   2 p accepted.
  | Species Mg :    hydro orbital   2 p accepted.
  | Species Mg :    ionic orbital   3 d accepted.
  | Species O :   atomic orbital   1 s accepted.
  | Species O :    hydro orbital   3 s accepted.
  | Species O :   atomic orbital   2 s accepted.
  | Species O :   atomic orbital   2 p accepted.
  | Species O :    hydro orbital   2 p accepted.
  | Species O :    hydro orbital   3 d accepted.

  Basis size parameters after reduction:
  | Total number of radial functions:       13
  | Total number of basis functions :     3132

  Per-task memory consumption for arrays in subroutine allocate_ext:
  |           4.515952MB.
  Testing on-site integration grid accuracy.
  |  Species  Function  <phi|h_atom|phi> (log., in eV)  <phi|h_atom|phi> (rad., in eV)
           1        1              -1261.4972086064              -1261.4966102663
           1        2                -79.6612686792                -79.6612661850
           1        3                 -4.5766449847                 -4.5764463354
           1        4                  4.2856803862                  4.2813226631
           1        5                -46.4234362291                -46.4234361907
           1        6                 -0.5282558506                 -0.5281575826
           1        7                  4.0144889854                  4.0133736295
           2        8               -515.0295626295               -515.0294562738
           2        9                 15.1698434419                 15.1699322204
           2       10                -21.6038822678                -21.6039118105
           2       11                 -9.0211123393                 -9.0212703513
           2       12                  8.3047696391                  8.2854840999
           2       13                 45.8428042125                 45.8427461222

  Preparing densities etc. for the partition functions (integrals / Hartree potential).

  Preparations completed.
  max(cpu_time)          :      0.217 s.
  Wall clock time (cpu1) :      0.415 s.
------------------------------------------------------------

************************** W A R N I N G *******************************
* Skipping the SCF initialization for now - done inside wrapper      *
************************************************************************

  @ DRIVER MODE: Connecting to host:port localhost       12345
  @ DRIVER MODE: Message from server: STATUS
  @ DRIVER MODE: Message from server: POSDATA
  @ DRIVER MODE: Received positions

------------------------------------------------------------
          Begin self-consistency loop: Initialization.

          Date     :  20191125, Time     :  161253.166
------------------------------------------------------------

  Initializing index lists of integration centers etc. from given atomic structure:
  Mapping all atomic coordinates to central unit cell.

  Initializing the k-points
  Using symmetry for reducing the k-points
  | k-points reduced from:        8 to        8
  | Number of k-points                             :         8
  The eigenvectors in the calculations are REAL.
  | Number of basis functions in the Hamiltonian integrals :     14787
  | Number of basis functions in a single unit cell        :      3132
  | Number of centers in hartree potential         :      3756
  | Number of centers in hartree multipole         :      3096
  | Number of centers in electron density summation:      2296
  | Number of centers in basis integrals           :      2572
  | Number of centers in integrals                 :      1215
  | Number of centers in hamiltonian               :      2680
  | Consuming          3 KiB for k_phase.
  | Number of super-cells (origin) [n_cells]                     :         125
  | Number of super-cells (after PM_index) [n_cells]             :          28
  | Number of super-cells in hamiltonian [n_cells_in_hamiltonian]:          28
  | Size of matrix packed + index [n_hamiltonian_matrix_size] :     7376421
  | Estimated reciprocal-space cutoff momentum G_max:         3.25739508 bohr^-1 .
  | Reciprocal lattice points for long-range Hartree potential:    8120
  Using simple linear distribution as default method for k-point parallelism.
  * Using 192 tasks for Scalapack Eigenvalue solver.
  Detailed listing of tasks and assigned k-points:
  ** shortened for testing
  ** ...

  What follows are estimated values for band gap, HOMO, LUMO, etc.
  | They are estimated on a discrete k-point grid and not necessarily exact.
  | For converged numbers, create a DOS and/or band structure plot on a denser k-grid.

  Highest occupied state (VBM) at    -12.42256568 eV (relative to internal zero)
  | Occupation number:      2.00000000
  | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

  Lowest unoccupied state (CBM) at    -6.03116942 eV (relative to internal zero)
  | Occupation number:      0.00000000
  | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

  ESTIMATED overall HOMO-LUMO gap:      6.39139626 eV between HOMO at k-point 1 and LUMO at k-point 1
  | This appears to be a direct band gap.
  The gap value is above 0.2 eV. Unless you are using a very sparse k-point grid,
  this system is most likely an insulator or a semiconductor.
  Calculating total energy contributions from superposition of free atom densities.

  Total energy components:
  | Sum of eigenvalues            :      -16390.95968681 Ha     -446020.70636178 eV
  | XC energy correction          :       -2629.08068281 Ha      -71540.92534156 eV
  | XC potential correction       :        3424.46246742 Ha       93184.36490682 eV
  | Free-atom electrostatic energy:      -14136.06164921 Ha     -384661.80885229 eV
  | Hartree energy correction     :           0.00000000 Ha           0.00000000 eV
  | Entropy correction            :           0.00000000 Ha           0.00000000 eV
  | ---------------------------
  | Total energy                  :      -29731.63955141 Ha     -809039.07564882 eV
  | Total energy, T -> 0          :      -29731.63955141 Ha     -809039.07564882 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy        :      -29731.63955141 Ha     -809039.07564882 eV

  Derived energy quantities:
  | Kinetic energy                :       29811.78541900 Ha      811219.95566799 eV
  | Electrostatic energy          :      -56914.34428760 Ha    -1548718.10597524 eV
  | Energy correction for multipole
  | error in Hartree potential    :           0.00000000 Ha           0.00000000 eV
  | Sum of eigenvalues per atom                           :       -2064.91067760 eV
  | Total energy (T->0) per atom                          :       -3745.55127615 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy per atom                       :       -3745.55127615 eV
  Initialize hartree_potential_storage
  Max. number of atoms included in rho_multipole:          216

  End scf initialization - timings             :  max(cpu_time)    wall_clock(cpu1)
  | Time for scf. initialization               :        9.052 s           9.165 s
  | Boundary condition initialization          :        1.528 s           1.531 s
  | Integration                                :        2.692 s           2.699 s
  | Solution of K.-S. eqns.                    :        1.866 s           1.874 s
  | Grid partitioning                          :        0.648 s           0.631 s
  | Preloading free-atom quantities on grid    :        0.006 s           0.010 s
  | Free-atom superposition energy             :        0.279 s           0.280 s
  | Total energy evaluation                    :        0.003 s           0.010 s

  Partial memory accounting:
  | Current value for overall tracked memory usage:
  |   Minimum:      113.209 MB (on task  19)
  |   Maximum:      113.543 MB (on task   0)
  |   Average:      113.423 MB
  | Peak value for overall tracked memory usage:
  |   Minimum:      225.764 MB (on task  19 after allocating temp_ham_ovlp)
  |   Maximum:      226.098 MB (on task   0 after allocating temp_ham_ovlp)
  |   Average:      225.978 MB
  | Largest tracked array allocation so far:
  |   Minimum:       56.278 MB (overlap_matrix on task   0)
  |   Maximum:       56.278 MB (overlap_matrix on task   0)
  |   Average:       56.278 MB
  Note:  These values currently only include a subset of arrays which are explicitly tracked.
  The "true" memory usage will be greater.
------------------------------------------------------------
Convergence:    q app. |  density  | eigen (eV) | Etot (eV) | forces (eV/A) |       CPU time |     Clock time
  SCF    1 :  0.54E-11 |  0.28E+01 |   0.12E+04 |  0.12E+03 |             . |        7.567 s |        7.595 s
  SCF    2 :  0.68E-11 |  0.12E+01 |   0.46E+03 |  0.25E+02 |             . |        7.857 s |        7.857 s
  SCF    3 :  0.56E-11 |  0.68E+00 |  -0.44E+03 |  0.41E+01 |             . |        8.146 s |        8.145 s
  SCF    4 :  0.43E-11 |  0.96E-01 |  -0.57E+01 | -0.94E-01 |             . |        7.457 s |        7.457 s
  SCF    5 :  0.52E-11 |  0.22E+00 |  -0.43E+01 |  0.20E+00 |             . |        8.219 s |        8.219 s
  SCF    6 :  0.46E-11 |  0.85E-01 |  -0.33E+01 | -0.26E-01 |             . |        7.422 s |        7.422 s
  SCF    7 :  0.50E-11 |  0.11E+00 |  -0.87E+01 |  0.70E-01 |             . |        7.328 s |        7.327 s
  SCF    8 :  0.47E-11 |  0.15E-01 |  -0.11E+01 |  0.29E-03 |             . |        7.478 s |        7.478 s
  SCF    9 :  0.45E-11 |  0.19E-01 |  -0.14E+01 |  0.20E-02 |             . |        7.306 s |        7.306 s
  SCF   10 :  0.61E-11 |  0.13E-01 |  -0.32E+01 |  0.35E-02 |             . |        7.319 s |        7.319 s
  SCF   90 :  0.57E-11 |  0.91E-06 |   0.16E-05 |  0.18E-08 |             . |        7.199 s |        7.199 s
  SCF   91 :  0.58E-11 |  0.16E-05 |   0.27E-05 |  0.21E-08 |      0.17E-01 |       29.010 s |       29.016 s

  Total energy components:
  | Sum of eigenvalues            :      -16347.85189592 Ha     -444847.68368895 eV
  | XC energy correction          :       -2639.97288985 Ha      -71837.31737527 eV
  | XC potential correction       :        3438.67297329 Ha       93571.05244586 eV
  | Free-atom electrostatic energy:      -14136.06164921 Ha     -384661.80885229 eV
  | Hartree energy correction     :         -40.81786389 Ha       -1110.71058866 eV
  | Entropy correction            :           0.00000000 Ha           0.00000000 eV
  | ---------------------------
  | Total energy                  :      -29726.03132558 Ha     -808886.46805931 eV
  | Total energy, T -> 0          :      -29726.03132558 Ha     -808886.46805931 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy        :      -29726.03132558 Ha     -808886.46805931 eV

  Derived energy quantities:
  | Kinetic energy                :       29783.32759545 Ha      810445.57888928 eV
  | Electrostatic energy          :      -56869.38603118 Ha    -1547494.72957331 eV
  | Energy correction for multipole
  | error in Hartree potential    :           0.37332979 Ha          10.15882051 eV
  | Sum of eigenvalues per atom                           :       -2059.48001708 eV
  | Total energy (T->0) per atom                          :       -3744.84475953 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy per atom                       :       -3744.84475953 eV
  What follows are estimated values for band gap, HOMO, LUMO, etc.
  | They are estimated on a discrete k-point grid and not necessarily exact.
  | For converged numbers, create a DOS and/or band structure plot on a denser k-grid.

  Highest occupied state (VBM) at     -9.65531162 eV (relative to internal zero)
  | Occupation number:      2.00000000
  | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

  Lowest unoccupied state (CBM) at    -5.16745431 eV (relative to internal zero)
  | Occupation number:      0.00000000
  | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

  ESTIMATED overall HOMO-LUMO gap:      4.48785731 eV between HOMO at k-point 1 and LUMO at k-point 1
  | This appears to be a direct band gap.
  The gap value is above 0.2 eV. Unless you are using a very sparse k-point grid,
  this system is most likely an insulator or a semiconductor.

  Self-consistency cycle converged.

  Removing unitary transformations (pure translations, rotations) from forces on atoms.
  Atomic forces before filtering:
  | Net force on center of mass :   0.177145E-02  0.850493E-03  0.241416E-02 eV/A
  Atomic forces after filtering:
  | Net force on center of mass :  -0.179800E-15 -0.564488E-16 -0.433820E-16 eV/A

  Energy and forces in a compact form:
  | Total energy uncorrected      :         -0.808886468059309E+06 eV
  | Total energy corrected        :         -0.808886468059309E+06 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy        :         -0.808886468059309E+06 eV
  Total atomic forces (unitary forces cleaned) [eV/Ang]:
  |    1         -0.119012361951726E+00          0.275516380987299E-05         -0.317395561009859E-04
  |    2         -0.518759077946052E-02         -0.485759862222911E-04         -0.565082705552965E-04
  |    3         -0.494679373402557E-03         -0.647456495281857E-05         -0.140350415740300E-04
  |    4         -0.379629476838036E-03          0.551687518626788E-05          0.305282337568159E-04
  |    5         -0.496457938464057E-03         -0.292444699386120E-05         -0.102259003287427E-04
  |    6         -0.518879823870005E-02          0.427152085393204E-04          0.349220747889230E-04
  |    7          0.532237610724758E-02         -0.167095450449332E-04          0.167237857887971E-01
  |    8         -0.368412132532592E-03          0.808116481442259E-04          0.119209141820651E-02
  |    9          0.813170065600783E-04          0.101319661272773E-03          0.321891746986132E-04
  |   10          0.138420662577077E-03          0.115262865841199E-04         -0.124679725703046E-03
  |   11         -0.369755936081618E-03         -0.121402052319909E-02         -0.955714043896540E-04
  |   12          0.531092903788260E-02         -0.167435519287347E-01         -0.254387866330507E-04
  |   13          0.470050033081477E-03         -0.327548611858284E-04          0.121228448366338E-02
  |   14          0.720105874057748E-03          0.895680127329521E-04          0.291244187423601E-04
  |   15          0.418577122626110E-03          0.366549581224780E-03         -0.354107870541968E-03
  |   16          0.687647911630069E-03          0.192945778087496E-04         -0.146258666971794E-03
  |   17          0.469074114115072E-03         -0.122188614310707E-02         -0.131855652937891E-04
  |   18          0.164131577646961E-02         -0.121336431152045E-02          0.118927005973412E-02
  |   19         -0.218990106337357E-03         -0.291221178460380E-04          0.221365175760942E-04
  |   20          0.138548532840371E-03          0.121612778852625E-04         -0.182495469915569E-04
  |   21          0.207329283052566E-03         -0.664845521831192E-05         -0.286395412960525E-04
  |   22         -0.109109670577258E-03          0.275057799822283E-04         -0.190074886382147E-05
  |   23          0.204672350161883E-03         -0.541800698187812E-05         -0.502070008352067E-05
  |   24          0.129288053614765E-03         -0.143172363004285E-04         -0.168136058947289E-04
  |   25          0.474925794621975E-03         -0.114873920220358E-04          0.121250401773246E-02
  |   26          0.164593777070683E-02         -0.121287638661786E-02          0.119048016337496E-02
  |   27          0.473987220058519E-03         -0.122273466759395E-02         -0.183024353842434E-04
  |   28          0.693287042368994E-03          0.238053030113415E-04         -0.155537966043552E-03
  |   29          0.414926707766132E-03          0.361373797741739E-03         -0.355819970971974E-03
  |   30          0.717322437858304E-03          0.101184590294972E-03          0.235446942562945E-04
  |   31          0.531645571542387E-02          0.104233715841785E-05          0.166211869543140E-01
  |   32          0.530597036676941E-02         -0.166401923022732E-01         -0.281015796427600E-04
  |   33         -0.361345330684274E-03         -0.121885456593814E-02         -0.100357697598852E-03
  |   34          0.131384926559678E-03          0.130156982510534E-04         -0.114511654923721E-03
  |   35          0.134391351616879E-03          0.958607110939483E-04          0.276877214814083E-04
  |   36         -0.361117112362966E-03          0.871032009116820E-04          0.119726875867392E-02
  |   37          0.531261071504286E-02          0.167367249299646E-01         -0.314077577575146E-04
  |   38         -0.368276448660070E-03          0.120588048346307E-02          0.731353427387136E-04
  |   39          0.140091617009063E-03          0.167838101908874E-04          0.857944193523855E-04
  |   40          0.943537929732435E-04         -0.110285343730211E-03          0.278358151803988E-04
  |   41         -0.370175598544466E-03         -0.876639366139021E-04         -0.121560534891810E-02
  |   42          0.532249110104152E-02         -0.115527626703236E-04         -0.167453821364556E-01
  |   43          0.164222762543980E-02          0.120536788418716E-02          0.118929198798053E-02
  |   44          0.419127262306640E-03          0.366969329491016E-03          0.328192056381104E-03
  |   45         -0.567429334360479E-03          0.175920146928126E-04          0.342031119116756E-04
  |   46          0.419253955119386E-03         -0.376264868332112E-03         -0.354547256935620E-03
  |   47          0.164170941878801E-02         -0.121283727407694E-02         -0.121091398921010E-02
  |   48          0.878218756550891E-02         -0.285368215070505E-04         -0.235010490292477E-04
  |   49          0.139922978190709E-03          0.112910775211164E-04         -0.146476724329435E-04
  |   50         -0.514882333400361E-03         -0.302538122768968E-04          0.429725485237623E-04
  |   51         -0.414465564108598E-03          0.177110557136156E-04         -0.134096251098733E-04
  |   52          0.204608474268078E-03         -0.356907815246535E-05         -0.308452939661014E-04
  |   53          0.251723191232588E-02         -0.372142786168174E-05         -0.131503690737709E-04
  |   54          0.252912993515159E-02         -0.324833166677647E-04         -0.107800109834371E-04
  |   55          0.722257568344988E-03         -0.105189440185260E-03          0.273347163261817E-04
  |   56          0.419715007470126E-03         -0.374369434061566E-03          0.326731754521336E-03
  |   57          0.696265497987473E-03          0.176862053592593E-04          0.115082238456575E-03
  |   58          0.475466724529142E-03          0.121453995422453E-02         -0.122426144366561E-04
  |   59          0.164780676956695E-02          0.120520286719222E-02         -0.121151080174965E-02
  |   60          0.475019973269784E-03         -0.152914184439519E-04         -0.123677472795820E-02
  |   61         -0.362126226433504E-03         -0.923275928333968E-04          0.119672953035639E-02
  |   62         -0.362156914169821E-03         -0.121904706534452E-02          0.789289495223576E-04
  |   63         -0.131837292674444E-03          0.165967815761392E-04         -0.232518014213093E-04
  |   64         -0.360785795295535E-03          0.120985864033755E-02         -0.999892126142993E-04
  |   65         -0.360741979017746E-03          0.868466855556252E-04         -0.122108750425235E-02
  |   66         -0.132774088803786E-03         -0.792500306523515E-05          0.277998512218076E-04
  |   67         -0.518836650813733E-02         -0.490348556922105E-04          0.349748865281067E-04
  |   68         -0.454338912837790E-03         -0.210726999938088E-04         -0.314366401906128E-04
  |   69         -0.341441961819019E-03          0.147596233698625E-04         -0.142625625497577E-04
  |   70         -0.496190583127518E-03         -0.337989499509331E-05         -0.135164013280220E-04
  |   71         -0.344748714219398E-03         -0.223476379688594E-05          0.294121538195963E-04
  |   72         -0.452286795590121E-03         -0.600968160442960E-05          0.500300654232741E-05
  |   73          0.470173952816009E-03          0.121371215399748E-02         -0.213225355994191E-04
  |   74          0.690481235882694E-03          0.270149054641106E-04          0.106771446550675E-03
  |   75          0.420227811172740E-03         -0.376940418491807E-03          0.328960088716528E-03
  |   76          0.723922265170705E-03         -0.924820308188315E-04          0.233835336582036E-04
  |   77          0.470607956261532E-03         -0.264742376640130E-04         -0.123618210608853E-02
  |   78          0.164297444688402E-02          0.120466914943541E-02         -0.121088713184674E-02
  |   79          0.208605640480606E-03         -0.583581176739221E-05         -0.463779205963306E-05
  |   80         -0.413246048867819E-03          0.249595982893479E-04         -0.535131724688666E-04
  |   81         -0.485064260207140E-03          0.929407484336317E-05          0.371349217274374E-04
  |   82          0.129285268393642E-03         -0.155516587945251E-04         -0.160974424504926E-04
  |   83          0.252964906260363E-02         -0.255012846067926E-04         -0.104215183353084E-04
  |   84          0.251658658671213E-02         -0.417268505267791E-05         -0.208859964236486E-04
  |   85          0.419326021915699E-03         -0.374620935875798E-03         -0.355288927180991E-03
  |   86         -0.592957231107338E-03          0.150857373311209E-04          0.332482515924373E-04
  |   87          0.415332412342989E-03          0.362562840103004E-03          0.327788892730229E-03
  |   88          0.164762691864482E-02          0.120534850337802E-02          0.119011678515242E-02
  |   89          0.877228721516628E-02         -0.147597853444878E-04         -0.214855183264954E-04
  |   90          0.164596436852324E-02         -0.121224162896196E-02         -0.121151604212109E-02
  |   91          0.843956985949821E-04         -0.109489293386209E-03          0.304769739451541E-04
  |   92          0.133969685429567E-03          0.914764182733951E-05          0.734929296262536E-04
  |   93         -0.361246633650066E-03          0.121030431858322E-02          0.783921423831232E-04
  |   94          0.530698510916603E-02          0.166338822427049E-01         -0.245707568383541E-04
  |   95          0.531708403241977E-02         -0.156506154915326E-05         -0.166423716389587E-01
  |   96         -0.361362706150616E-03         -0.917746510749925E-04         -0.122086190547726E-02
  |   97         -0.495467311507061E-03         -0.692508998836333E-05         -0.982880329504751E-05
  |   98         -0.341987098680571E-03          0.104110483987335E-04         -0.117776084491391E-04
  |   99         -0.454817935986620E-03          0.128311197216138E-04         -0.271311671636539E-04
  |  100         -0.518855582713648E-02          0.424275057441895E-04         -0.564911570852096E-04
  |  101         -0.451793359084410E-03         -0.790103691787807E-05         -0.290258078336556E-04
  |  102         -0.343997349716902E-03         -0.538331099530636E-05          0.312974881942226E-04
  |  103         -0.368499498946202E-03          0.120547036337704E-02         -0.953010474623487E-04
  |  104         -0.136203072998444E-03          0.175785385856787E-04         -0.250910214110154E-04
  |  105         -0.369772143640060E-03         -0.121443019558008E-02          0.735238147570348E-04
  |  106         -0.370314462983028E-03         -0.881534018508957E-04          0.119206783005834E-02
  |  107         -0.136873132793777E-03         -0.218118906445721E-04          0.281216449362154E-04
  |  108         -0.368446644653551E-03          0.802989220290471E-04         -0.121563370434308E-02
  |  109          0.455479334134689E-03         -0.318754151368766E-02         -0.319466901655521E-02
  |  110          0.166138264605238E-03         -0.465094809008625E-03         -0.468986339576328E-03
  |  111          0.686757419883299E-04          0.723433231428553E-04          0.203598131431912E-03
  |  112          0.164832726693664E-03          0.457096441618671E-03          0.445756463636897E-03
  |  113          0.454647579018006E-03          0.317994105464319E-02          0.317163065613774E-02
  |  114          0.161205169749306E-01         -0.493995121088980E-04         -0.116117353283496E-03
  |  115         -0.522572425436770E-03         -0.288817404236280E-03         -0.630116126175970E-03
  |  116         -0.122541028821937E-03         -0.180699973833249E-03          0.212789143112284E-03
  |  117         -0.896833428875549E-04          0.994905150487619E-04          0.143825233858540E-03
  |  118         -0.520999369338172E-03          0.619771047825692E-03          0.272048390350830E-03
  |  119         -0.335812219869158E-02          0.289395497462209E-02         -0.511050725292329E-04
  |  120         -0.335816246784530E-02         -0.162525320264281E-03         -0.290853821344318E-02
  |  121         -0.714176864298821E-03          0.172109906670655E-03          0.186177561035286E-03
  |  122         -0.180073125931452E-03         -0.341764278761372E-04         -0.785020867687771E-04
  |  123         -0.781413640611539E-03          0.132797803381638E-03         -0.198328955518428E-03
  |  124         -0.103538682954804E-02         -0.381296252477668E-05          0.320647021149128E-05
  |  125         -0.152832786117954E-02         -0.445104950581878E-05         -0.110567679836288E-04
  |  126         -0.103593941235340E-02         -0.146434577471333E-03         -0.116388916999341E-04
  |  127         -0.526104252116470E-03          0.281561718400453E-03         -0.630423754331310E-03
  |  128         -0.524368947200821E-03          0.619367978288964E-03         -0.295672516029959E-03
  |  129         -0.212204733457306E-03          0.136721846123003E-03         -0.269420390288296E-04
  |  130         -0.523647019091070E-03         -0.627780437565774E-03          0.273553692373832E-03
  |  131         -0.524459187459343E-03         -0.289633249659508E-03          0.607979966170910E-03
  |  132         -0.223973520598026E-03         -0.479694362161180E-04          0.156207731708732E-03
  |  133          0.429204361517739E-03          0.317691846897217E-02         -0.319190485670363E-02
  |  134          0.803150558740872E-03          0.753283575354350E-03         -0.988931754242789E-04
  |  135          0.185046398781922E-03          0.100440417066264E-03          0.975511939529324E-04
  |  136          0.161161883947988E-03         -0.466732571411339E-03          0.447918405653638E-03
  |  137          0.182530432944151E-03         -0.112771709092091E-03          0.179198621850356E-03
  |  138          0.803468985712579E-03         -0.414752200814913E-05         -0.767417738963475E-03
  |  139          0.143206881679199E-01         -0.330710283240411E-04         -0.131111031085756E-03
  |  140          0.999244135634077E-03          0.389409086824781E-05         -0.452468832766340E-05
  |  141          0.313916281933964E-03          0.677391227560589E-04         -0.120599057679184E-04
  |  142          0.318063506847095E-03         -0.405832076664253E-05          0.189601965477551E-03
  |  143          0.996995205073439E-03         -0.105285519711170E-04         -0.193473808028341E-04
  |  144          0.143206463601523E-01          0.235897890691540E-04          0.176061641000115E-04
  |  145         -0.520060241568279E-03         -0.627178642589187E-03         -0.294890908354749E-03
  |  146         -0.879266979640786E-04          0.137724060031698E-03         -0.194115341887655E-03
  |  147         -0.117392211861083E-03          0.170983789997879E-03          0.176261837011870E-03
  |  148         -0.522909317281839E-03          0.281008935146670E-03          0.607390077920003E-03
  |  149         -0.335786877137545E-02         -0.118995670705305E-03          0.288577550134907E-02
  |  150         -0.335767398116911E-02         -0.290206367320085E-02         -0.974937638502049E-04
  |  151         -0.177573401470858E-03         -0.397387880492195E-04         -0.121306493262030E-04
  |  152         -0.133370516247777E-03          0.997741225433267E-04          0.221402352768127E-03
  |  153         -0.191503280982715E-03          0.301133393126951E-04         -0.579861286240614E-04
  |  154         -0.152830573574579E-02         -0.405922815767363E-05         -0.115911073746599E-04
  |  155         -0.281024711034494E-02         -0.153702896654514E-03         -0.636915749409922E-04
  |  156         -0.152885185551779E-02         -0.350239516656068E-05         -0.109846410176659E-04
  |  157         -0.135697797154595E-03          0.189710702210611E-03          0.200887202520731E-03
  |  158         -0.824766153237564E-04          0.747678546829478E-04         -0.229789676852593E-03
  |  159         -0.523544791449824E-03         -0.627409904300937E-03         -0.296150391364797E-03
  |  160         -0.336046484086929E-02         -0.290573986988598E-02         -0.471153096334804E-04
  |  161         -0.336092198554407E-02         -0.771440674672303E-04          0.288994460547350E-02
  |  162         -0.526011511200246E-03          0.281804753165727E-03          0.607529226428980E-03
  |  163          0.161856524297386E-03          0.457465246192158E-03         -0.469727421312690E-03
  |  164          0.184822539043206E-03          0.774671694926839E-04         -0.120241589687033E-03
  |  165          0.803605968266618E-03         -0.761001675838423E-03         -0.767848028439445E-04
  |  166          0.429897648473271E-03         -0.318498407253961E-02          0.316943262248778E-02
  |  167          0.803751981823900E-03         -0.118962882522971E-04          0.744890832752938E-03
  |  168          0.183126522528060E-03          0.104524084284549E-03          0.187613468090597E-03
  |  169          0.998607371973313E-03          0.364458065665340E-05         -0.180559256307077E-04
  |  170          0.550902997925187E-03          0.894926504484323E-04         -0.115266302057690E-03
  |  171          0.998587765573283E-03         -0.118258186153997E-04         -0.422772792376669E-05
  |  172          0.997125703531043E-03         -0.107307783580050E-04         -0.301258495015007E-05
  |  173          0.553610625282384E-03         -0.416093005434179E-04          0.197468485453838E-03
  |  174          0.997361982843538E-03          0.274604423743994E-05         -0.189976352745187E-04
  |  175          0.819233266351670E-03         -0.750665683610301E-03         -0.119524229657576E-03
  |  176          0.189936512306804E-03          0.119731652693511E-03         -0.117230101942459E-03
  |  177          0.164934709649158E-03          0.456633146257920E-03         -0.468193103452903E-03
  |  178          0.187302829617291E-03          0.101395172425963E-03          0.177630789607491E-03
  |  179          0.819049160507982E-03         -0.732278868328353E-04          0.734131503201529E-03
  |  180          0.455409192417670E-03         -0.318795643676731E-02          0.317176351572238E-02
  |  181         -0.768503302746111E-03          0.175072098601091E-03          0.607632777894234E-04
  |  182         -0.191832944786343E-03          0.416042804907856E-04         -0.421006225136861E-04
  |  183         -0.782934542048400E-03         -0.159048630567893E-03          0.153218774913745E-03
  |  184         -0.103559585776758E-02         -0.111600652334010E-03         -0.115315120724235E-04
  |  185         -0.152840085369083E-02         -0.401929596045748E-05         -0.116393126644795E-04
  |  186         -0.103563840328247E-02         -0.442215216060914E-05         -0.406796719830890E-04
  |  187         -0.837661739816599E-04          0.101378335567555E-03          0.173089045650909E-03
  |  188         -0.211479991719774E-03         -0.207278254906739E-03          0.176353634003782E-03
  |  189         -0.524474397198791E-03         -0.289512985205629E-03         -0.630871288651830E-03
  |  190         -0.336090285157812E-02         -0.533099518927347E-04         -0.291263947847707E-02
  |  191         -0.336118217433320E-02          0.289785010396676E-02         -0.740745963450276E-04
  |  192         -0.524513837598162E-03          0.619677626477388E-03          0.273166744252878E-03
  |  193         -0.659575919510912E-05          0.414713449738465E-04          0.192804600515189E-03
  |  194          0.161034923710800E-03         -0.466434550715283E-03         -0.470542886951697E-03
  |  195          0.429841441558390E-03         -0.318472955985017E-02         -0.319221386838936E-02
  |  196          0.139887835512905E-01          0.227682118537222E-04         -0.997032006522022E-04
  |  197          0.429679746998123E-03          0.317721401711178E-02          0.316911332113803E-02
  |  198          0.162388931719975E-03          0.457720914497143E-03          0.447059582061974E-03
  |  199          0.313508450038135E-03          0.465013898173586E-04         -0.140572008604119E-04
  |  200          0.998397092059462E-03         -0.114884836351944E-04         -0.184963505291876E-04
  |  201          0.143204221045434E-01          0.253994270188708E-04         -0.109661751934460E-03
  |  202          0.143209682221198E-01          0.549185112191198E-05         -0.405424357484979E-04
  |  203          0.997902882562226E-03          0.311206287177173E-05         -0.347123971448758E-05
  |  204          0.330957484551552E-03         -0.436241586930858E-05          0.207190041447451E-03
  |  205          0.189651639940363E-03          0.836995338461696E-04          0.944607863481205E-04
  |  206          0.818481535034793E-03          0.742975512966809E-03         -0.826306048576369E-04
  |  207          0.454721148564881E-03          0.317965592023785E-02         -0.319459179210139E-02
  |  208          0.819134594315545E-03         -0.987249062021017E-04         -0.756572551019718E-03
  |  209          0.188390047059969E-03         -0.109534843824609E-03          0.202040404477769E-03
  |  210          0.166008326906956E-03         -0.465473859931084E-03          0.446344933293556E-03
  |  211         -0.209435689781824E-03          0.142245052408247E-03         -0.451706684219384E-04
  |  212         -0.520369965919543E-03          0.619141991462990E-03         -0.294781783467021E-03
  |  213         -0.523050250367518E-03          0.280921584353888E-03         -0.630194161838450E-03
  |  214         -0.221190054151161E-03         -0.155053319546961E-03          0.160993726328277E-03
  |  215         -0.522016570308690E-03         -0.289045415328539E-03          0.607263245992090E-03
  |  216         -0.520288781674739E-03         -0.627733279543286E-03          0.272478097279923E-03

  Removing unitary transformations (pure translations, rotations) from forces on atoms.
  Atomic forces before filtering:
  | Net force on center of mass :  -0.179800E-15 -0.564488E-16 -0.433820E-16 eV/A
  Atomic forces after filtering:
  | Net force on center of mass :  -0.179626E-15 -0.831923E-16 -0.429900E-16 eV/A
  @ DRIVER - WARNING: the stress tensor will not be computed
 ------------------------------------------------------------------------
 Atomic structure that was used in the preceding time step of the wrapper
                         x [A]             y [A]             z [A]
  lattice_vector        12.67881443        0.00000000        0.00000000
  lattice_vector         0.00000000       12.67881443        0.00000000
  lattice_vector        -0.00000000       -0.00000000       12.67881443

            atom         0.01000000       -0.00000000        0.00000000  Mg
            atom        -0.00000000        2.11313574        2.11313574  Mg
            atom        -0.00000000        4.22627148        4.22627148  Mg
            atom        -0.00000000        6.33940721        6.33940721  Mg
            atom        -0.00000000       -4.22627148       -4.22627148  Mg
            atom        -0.00000000       -2.11313574       -2.11313574  Mg
            atom         2.11313574       -0.00000000        2.11313574  Mg
            atom         2.11313574        2.11313574        4.22627148  Mg
            atom         2.11313574        4.22627148        6.33940721  Mg
            atom         2.11313574        6.33940721       -4.22627148  Mg
            atom         2.11313574       -4.22627148       -2.11313574  Mg
            atom         2.11313574       -2.11313574       -0.00000000  Mg
            atom         4.22627148       -0.00000000        4.22627148  Mg
            atom         4.22627148        2.11313574        6.33940721  Mg
            atom         4.22627148        4.22627148       -4.22627148  Mg
            atom         4.22627148        6.33940721       -2.11313574  Mg
            atom         4.22627148       -4.22627148       -0.00000000  Mg
            atom         4.22627148       -2.11313574        2.11313574  Mg
            atom         6.33940721       -0.00000000        6.33940721  Mg
            atom         6.33940721        2.11313574       -4.22627148  Mg
            atom         6.33940721        4.22627148       -2.11313574  Mg
            atom         6.33940721        6.33940721       -0.00000000  Mg
            atom         6.33940721       -4.22627148        2.11313574  Mg
            atom         6.33940721       -2.11313574        4.22627148  Mg
            atom        -4.22627148       -0.00000000       -4.22627148  Mg
            atom        -4.22627148        2.11313574       -2.11313574  Mg
            atom        -4.22627148        4.22627148       -0.00000000  Mg
            atom        -4.22627148        6.33940721        2.11313574  Mg
            atom        -4.22627148       -4.22627148        4.22627148  Mg
            atom        -4.22627148       -2.11313574        6.33940721  Mg
            atom        -2.11313574       -0.00000000       -2.11313574  Mg
            atom        -2.11313574        2.11313574        0.00000000  Mg
            atom        -2.11313574        4.22627148        2.11313574  Mg
            atom        -2.11313574        6.33940721        4.22627148  Mg
            atom        -2.11313574       -4.22627148        6.33940721  Mg
            atom        -2.11313574       -2.11313574       -4.22627148  Mg
            atom         2.11313574        2.11313574       -0.00000000  Mg
            atom         2.11313574        4.22627148        2.11313574  Mg
            atom         2.11313574        6.33940721        4.22627148  Mg
            atom         2.11313574       -4.22627148        6.33940721  Mg
            atom         2.11313574       -2.11313574       -4.22627148  Mg
            atom         2.11313574       -0.00000000       -2.11313574  Mg
            atom         4.22627148        2.11313574        2.11313574  Mg
            atom         4.22627148        4.22627148        4.22627148  Mg
            atom         4.22627148        6.33940721        6.33940721  Mg
            atom         4.22627148       -4.22627148       -4.22627148  Mg
            atom         4.22627148       -2.11313574       -2.11313574  Mg
            atom         4.22627148       -0.00000000        0.00000000  Mg
            atom         6.33940721        2.11313574        4.22627148  Mg
            atom         6.33940721        4.22627148        6.33940721  Mg
            atom         6.33940721        6.33940721       -4.22627148  Mg
            atom         6.33940721       -4.22627148       -2.11313574  Mg
            atom         6.33940721       -2.11313574       -0.00000000  Mg
            atom         6.33940721       -0.00000000        2.11313574  Mg
            atom        -4.22627148        2.11313574        6.33940721  Mg
            atom        -4.22627148        4.22627148       -4.22627148  Mg
            atom        -4.22627148        6.33940721       -2.11313574  Mg
            atom        -4.22627148       -4.22627148       -0.00000000  Mg
            atom        -4.22627148       -2.11313574        2.11313574  Mg
            atom        -4.22627148       -0.00000000        4.22627148  Mg
            atom        -2.11313574        2.11313574       -4.22627148  Mg
            atom        -2.11313574        4.22627148       -2.11313574  Mg
            atom        -2.11313574        6.33940721       -0.00000000  Mg
            atom        -2.11313574       -4.22627148        2.11313574  Mg
            atom        -2.11313574       -2.11313574        4.22627148  Mg
            atom        -2.11313574       -0.00000000        6.33940721  Mg
            atom         0.00000000        2.11313574       -2.11313574  Mg
            atom        -0.00000000        4.22627148       -0.00000000  Mg
            atom        -0.00000000        6.33940721        2.11313574  Mg
            atom        -0.00000000       -4.22627148        4.22627148  Mg
            atom        -0.00000000       -2.11313574        6.33940721  Mg
            atom         0.00000000       -0.00000000       -4.22627148  Mg
            atom         4.22627148        4.22627148       -0.00000000  Mg
            atom         4.22627148        6.33940721        2.11313574  Mg
            atom         4.22627148       -4.22627148        4.22627148  Mg
            atom         4.22627148       -2.11313574        6.33940721  Mg
            atom         4.22627148       -0.00000000       -4.22627148  Mg
            atom         4.22627148        2.11313574       -2.11313574  Mg
            atom         6.33940721        4.22627148        2.11313574  Mg
            atom         6.33940721        6.33940721        4.22627148  Mg
            atom         6.33940721       -4.22627148        6.33940721  Mg
            atom         6.33940721       -2.11313574       -4.22627148  Mg
            atom         6.33940721       -0.00000000       -2.11313574  Mg
            atom         6.33940721        2.11313574       -0.00000000  Mg
            atom        -4.22627148        4.22627148        4.22627148  Mg
            atom        -4.22627148        6.33940721        6.33940721  Mg
            atom        -4.22627148       -4.22627148       -4.22627148  Mg
            atom        -4.22627148       -2.11313574       -2.11313574  Mg
            atom        -4.22627148       -0.00000000        0.00000000  Mg
            atom        -4.22627148        2.11313574        2.11313574  Mg
            atom        -2.11313574        4.22627148        6.33940721  Mg
            atom        -2.11313574        6.33940721       -4.22627148  Mg
            atom        -2.11313574       -4.22627148       -2.11313574  Mg
            atom        -2.11313574       -2.11313574       -0.00000000  Mg
            atom        -2.11313574       -0.00000000        2.11313574  Mg
            atom        -2.11313574        2.11313574        4.22627148  Mg
            atom        -0.00000000        4.22627148       -4.22627148  Mg
            atom        -0.00000000        6.33940721       -2.11313574  Mg
            atom        -0.00000000       -4.22627148       -0.00000000  Mg
            atom        -0.00000000       -2.11313574        2.11313574  Mg
            atom         0.00000000       -0.00000000        4.22627148  Mg
            atom        -0.00000000        2.11313574        6.33940721  Mg
            atom         2.11313574        4.22627148       -2.11313574  Mg
            atom         2.11313574        6.33940721       -0.00000000  Mg
            atom         2.11313574       -4.22627148        2.11313574  Mg
            atom         2.11313574       -2.11313574        4.22627148  Mg
            atom         2.11313574       -0.00000000        6.33940721  Mg
            atom         2.11313574        2.11313574       -4.22627148  Mg
            atom         2.11313574        2.11313574        2.11313574  O
            atom         2.11313574        4.22627148        4.22627148  O
            atom         2.11313574        6.33940721        6.33940721  O
            atom         2.11313574       -4.22627148       -4.22627148  O
            atom         2.11313574       -2.11313574       -2.11313574  O
            atom         2.11313574        0.00000000       -0.00000000  O
            atom         4.22627148        2.11313574        4.22627148  O
            atom         4.22627148        4.22627148        6.33940721  O
            atom         4.22627148        6.33940721       -4.22627148  O
            atom         4.22627148       -4.22627148       -2.11313574  O
            atom         4.22627148       -2.11313574       -0.00000000  O
            atom         4.22627148        0.00000000        2.11313574  O
            atom         6.33940721        2.11313574        6.33940721  O
            atom         6.33940721        4.22627148       -4.22627148  O
            atom         6.33940721        6.33940721       -2.11313574  O
            atom         6.33940721       -4.22627148       -0.00000000  O
            atom         6.33940721       -2.11313574        2.11313574  O
            atom         6.33940721        0.00000000        4.22627148  O
            atom        -4.22627148        2.11313574       -4.22627148  O
            atom        -4.22627148        4.22627148       -2.11313574  O
            atom        -4.22627148        6.33940721       -0.00000000  O
            atom        -4.22627148       -4.22627148        2.11313574  O
            atom        -4.22627148       -2.11313574        4.22627148  O
            atom        -4.22627148        0.00000000        6.33940721  O
            atom        -2.11313574        2.11313574       -2.11313574  O
            atom        -2.11313574        4.22627148       -0.00000000  O
            atom        -2.11313574        6.33940721        2.11313574  O
            atom        -2.11313574       -4.22627148        4.22627148  O
            atom        -2.11313574       -2.11313574        6.33940721  O
            atom        -2.11313574        0.00000000       -4.22627148  O
            atom        -0.00000000        2.11313574       -0.00000000  O
            atom        -0.00000000        4.22627148        2.11313574  O
            atom        -0.00000000        6.33940721        4.22627148  O
            atom        -0.00000000       -4.22627148        6.33940721  O
            atom        -0.00000000       -2.11313574       -4.22627148  O
            atom        -0.00000000        0.00000000       -2.11313574  O
            atom         4.22627148        4.22627148        2.11313574  O
            atom         4.22627148        6.33940721        4.22627148  O
            atom         4.22627148       -4.22627148        6.33940721  O
            atom         4.22627148       -2.11313574       -4.22627148  O
            atom         4.22627148       -0.00000000       -2.11313574  O
            atom         4.22627148        2.11313574       -0.00000000  O
            atom         6.33940721        4.22627148        4.22627148  O
            atom         6.33940721        6.33940721        6.33940721  O
            atom         6.33940721       -4.22627148       -4.22627148  O
            atom         6.33940721       -2.11313574       -2.11313574  O
            atom         6.33940721       -0.00000000        0.00000000  O
            atom         6.33940721        2.11313574        2.11313574  O
            atom        -4.22627148        4.22627148        6.33940721  O
            atom        -4.22627148        6.33940721       -4.22627148  O
            atom        -4.22627148       -4.22627148       -2.11313574  O
            atom        -4.22627148       -2.11313574       -0.00000000  O
            atom        -4.22627148       -0.00000000        2.11313574  O
            atom        -4.22627148        2.11313574        4.22627148  O
            atom        -2.11313574        4.22627148       -4.22627148  O
            atom        -2.11313574        6.33940721       -2.11313574  O
            atom        -2.11313574       -4.22627148       -0.00000000  O
            atom        -2.11313574       -2.11313574        2.11313574  O
            atom        -2.11313574       -0.00000000        4.22627148  O
            atom        -2.11313574        2.11313574        6.33940721  O
            atom        -0.00000000        4.22627148       -2.11313574  O
            atom        -0.00000000        6.33940721       -0.00000000  O
            atom        -0.00000000       -4.22627148        2.11313574  O
            atom        -0.00000000       -2.11313574        4.22627148  O
            atom         0.00000000       -0.00000000        6.33940721  O
            atom         0.00000000        2.11313574       -4.22627148  O
            atom         2.11313574        4.22627148       -0.00000000  O
            atom         2.11313574        6.33940721        2.11313574  O
            atom         2.11313574       -4.22627148        4.22627148  O
            atom         2.11313574       -2.11313574        6.33940721  O
            atom         2.11313574       -0.00000000       -4.22627148  O
            atom         2.11313574        2.11313574       -2.11313574  O
            atom         6.33940721        6.33940721        2.11313574  O
            atom         6.33940721       -4.22627148        4.22627148  O
            atom         6.33940721       -2.11313574        6.33940721  O
            atom         6.33940721       -0.00000000       -4.22627148  O
            atom         6.33940721        2.11313574       -2.11313574  O
            atom         6.33940721        4.22627148       -0.00000000  O
            atom        -4.22627148        6.33940721        4.22627148  O
            atom        -4.22627148       -4.22627148        6.33940721  O
            atom        -4.22627148       -2.11313574       -4.22627148  O
            atom        -4.22627148       -0.00000000       -2.11313574  O
            atom        -4.22627148        2.11313574       -0.00000000  O
            atom        -4.22627148        4.22627148        2.11313574  O
            atom        -2.11313574        6.33940721        6.33940721  O
            atom        -2.11313574       -4.22627148       -4.22627148  O
            atom        -2.11313574       -2.11313574       -2.11313574  O
            atom        -2.11313574       -0.00000000        0.00000000  O
            atom        -2.11313574        2.11313574        2.11313574  O
            atom        -2.11313574        4.22627148        4.22627148  O
            atom        -0.00000000        6.33940721       -4.22627148  O
            atom        -0.00000000       -4.22627148       -2.11313574  O
            atom        -0.00000000       -2.11313574       -0.00000000  O
            atom         0.00000000       -0.00000000        2.11313574  O
            atom        -0.00000000        2.11313574        4.22627148  O
            atom        -0.00000000        4.22627148        6.33940721  O
            atom         2.11313574        6.33940721       -2.11313574  O
            atom         2.11313574       -4.22627148       -0.00000000  O
            atom         2.11313574       -2.11313574        2.11313574  O
            atom         2.11313574       -0.00000000        4.22627148  O
            atom         2.11313574        2.11313574        6.33940721  O
            atom         2.11313574        4.22627148       -4.22627148  O
            atom         4.22627148        6.33940721       -0.00000000  O
            atom         4.22627148       -4.22627148        2.11313574  O
            atom         4.22627148       -2.11313574        4.22627148  O
            atom         4.22627148       -0.00000000        6.33940721  O
            atom         4.22627148        2.11313574       -4.22627148  O
            atom         4.22627148        4.22627148       -2.11313574  O
 ------------------------------------------------------------------------
  @ DRIVER MODE: Message from server: STATUS
  @ DRIVER MODE: Message from server: GETFORCE
  @ DRIVER MODE: Returning v,forces,stress
  @ DRIVER MODE: Message from server: STATUS
  @ DRIVER MODE: Message from server: INIT
  @ DRIVER MODE: Receiving replica           0
  @ DRIVER MODE: Message from server: STATUS
  @ DRIVER MODE: Message from server: POSDATA
  @ DRIVER MODE: Received positions

------------------------------------------------------------
          Begin self-consistency loop: Re-initialization.

  Date     :  20191125, Time     :  162443.202
------------------------------------------------------------

  Initializing index lists of integration centers etc. from given atomic structure:
  Mapping all atomic coordinates to central unit cell.

  Initializing the k-points
  Using symmetry for reducing the k-points
  | k-points reduced from:        8 to        8
  | Number of k-points                             :         8
  The eigenvectors in the calculations are REAL.
  | Number of basis functions in the Hamiltonian integrals :     14787
  | Number of basis functions in a single unit cell        :      3132
  | Number of centers in hartree potential         :      3756
  | Number of centers in hartree multipole         :      3096
  | Number of centers in electron density summation:      2296
  | Number of centers in basis integrals           :      2572
  | Number of centers in integrals                 :      1215
  | Number of centers in hamiltonian               :      2680
  | Consuming          3 KiB for k_phase.
  | Number of super-cells (origin) [n_cells]                     :         125
  | Number of super-cells (after PM_index) [n_cells]             :          28
  | Number of super-cells in hamiltonian [n_cells_in_hamiltonian]:          28
  | Size of matrix packed + index [n_hamiltonian_matrix_size] :     7376421
  Partitioning the integration grid into batches with parallel hashing+maxmin method.
  Initializing partition tables, free-atom densities, potentials, etc. across the integration grid (initialize_grid_storage).
  | Species        1: outer_partition_radius set to              5.555717568450569 AA .
  | Species        2: outer_partition_radius set to              5.048384829883283 AA .
  | The sparse table of interatomic distances needs      12359.44 kbyte instead of     52921.47 kbyte of memory.
  | Net number of integration points:  1289088
  | of which are non-zero points    :   990776
  Renormalizing the initial density to the exact electron count on the 3D integration grid.
  | Initial density: Formal number of electrons (from input files) :    2160.0000000000
  | Integrated number of electrons on 3D grid     :    2160.3575811185
  | Charge integration error                      :       0.3575811185
  | Normalization factor for density and gradient :       0.9998344806
  Renormalizing the free-atom superposition density to the exact electron count on the 3D integration grid.
  | Formal number of electrons (from input files) :    2160.0000000000
  | Integrated number of electrons on 3D grid     :    2160.3575811185
  | Charge integration error                      :       0.3575811185
  | Normalization factor for density and gradient :       0.9998344806
  Calculating total energy contributions from superposition of free atom densities.
  Initialize hartree_potential_storage
  Max. number of atoms included in rho_multipole:          216
  Integrating overlap matrix.
  Time summed over all CPUs for integration: real work      178.859 s, elapsed      186.028 s
  Normalizing  ScaLAPACK eigenvectors

  End scf reinitialization - timings           :  max(cpu_time)    wall_clock(cpu1)
  | Time for scf. reinitialization             :        5.633 s           5.636 s
  | Boundary condition initialization          :        1.485 s           1.485 s
  | Integration                                :        1.401 s           1.356 s
  | Grid partitioning                          :        0.550 s           0.547 s
  | Preloading free-atom quantities on grid    :        1.540 s           1.539 s
  | Free-atom superposition energy             :        0.182 s           0.183 s
  | K.-S. eigenvector reorthonormalization     :        0.470 s           0.445 s
------------------------------------------------------------
Convergence:    q app. |  density  | eigen (eV) | Etot (eV) | forces (eV/A) |       CPU time |     Clock time
  SCF    1 :  0.48E-11 |  0.13E+01 |  -0.44E+06 | -0.81E+06 |             . |        6.972 s |        6.972 s
  SCF    2 :  0.36E-11 |  0.20E-02 |  -0.12E-02 |  0.13E-03 |             . |        7.247 s |        7.247 s
  SCF    3 :  0.43E-11 |  0.88E-03 |  -0.12E-02 |  0.29E-04 |             . |        7.239 s |        7.240 s
  SCF    4 :  0.42E-11 |  0.36E-03 |  -0.11E-02 |  0.15E-05 |             . |        7.242 s |        7.242 s
  SCF    5 :  0.40E-11 |  0.83E-04 |   0.37E-03 | -0.20E-06 |             . |        7.274 s |        7.274 s
  SCF    6 :  0.41E-11 |  0.29E-04 |   0.15E-03 | -0.26E-07 |             . |        7.240 s |        7.240 s
  SCF    7 :  0.41E-11 |  0.12E-04 |   0.82E-04 |  0.13E-07 |             . |        7.220 s |        7.220 s
  SCF    8 :  0.35E-11 |  0.36E-05 |   0.36E-04 |  0.99E-08 |             . |        7.241 s |        7.240 s
  SCF    9 :  0.45E-11 |  0.11E-05 |   0.91E-05 |  0.99E-09 |             . |        7.256 s |        7.256 s
  SCF   10 :  0.40E-11 |  0.41E-06 |   0.17E-05 |  0.48E-08 |             . |        7.209 s |        7.209 s
  SCF   11 :  0.46E-11 |  0.12E-06 |  -0.27E-06 | -0.79E-09 |      0.22E-01 |       29.143 s |       29.144 s

  Total energy components:
  | Sum of eigenvalues            :      -16347.85173669 Ha     -444847.67935619 eV
  | XC energy correction          :       -2639.97288881 Ha      -71837.31734700 eV
  | XC potential correction       :        3438.67297188 Ha       93571.05240749 eV
  | Free-atom electrostatic energy:      -14136.06170596 Ha     -384661.81039660 eV
  | Hartree energy correction     :         -40.81796663 Ha       -1110.71338460 eV
  | Entropy correction            :           0.00000000 Ha           0.00000000 eV
  | ---------------------------
  | Total energy                  :      -29726.03132622 Ha     -808886.46807690 eV
  | Total energy, T -> 0          :      -29726.03132622 Ha     -808886.46807690 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy        :      -29726.03132622 Ha     -808886.46807690 eV

  Derived energy quantities:
  | Kinetic energy                :       29783.32759108 Ha      810445.57877035 eV
  | Electrostatic energy          :      -56869.38602849 Ha    -1547494.72950025 eV
  | Energy correction for multipole
  | error in Hartree potential    :           0.37332967 Ha          10.15881726 eV
  | Sum of eigenvalues per atom                           :       -2059.47999702 eV
  | Total energy (T->0) per atom                          :       -3744.84475962 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy per atom                       :       -3744.84475962 eV
  What follows are estimated values for band gap, HOMO, LUMO, etc.
  | They are estimated on a discrete k-point grid and not necessarily exact.
  | For converged numbers, create a DOS and/or band structure plot on a denser k-grid.

  Highest occupied state (VBM) at     -9.65531857 eV (relative to internal zero)
  | Occupation number:      2.00000000
  | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

  Lowest unoccupied state (CBM) at    -5.16746014 eV (relative to internal zero)
  | Occupation number:      0.00000000
  | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

  ESTIMATED overall HOMO-LUMO gap:      4.48785843 eV between HOMO at k-point 1 and LUMO at k-point 1
  | This appears to be a direct band gap.
  The gap value is above 0.2 eV. Unless you are using a very sparse k-point grid,
  this system is most likely an insulator or a semiconductor.

  Self-consistency cycle converged.

  Removing unitary transformations (pure translations, rotations) from forces on atoms.
  Atomic forces before filtering:
  | Net force on center of mass :  -0.139592E-02  0.469618E-07 -0.602965E-07 eV/A
  Atomic forces after filtering:
  | Net force on center of mass :  -0.620240E-16  0.172482E-16 -0.682961E-16 eV/A

  Energy and forces in a compact form:
  | Total energy uncorrected      :         -0.808886468076903E+06 eV
  | Total energy corrected        :         -0.808886468076903E+06 eV  <-- do not rely on this value for anything but (periodic) metals
  | Electronic free energy        :         -0.808886468076903E+06 eV
  Total atomic forces (unitary forces cleaned) [eV/Ang]:
  |    1          0.454895003232345E-03         -0.318542358999639E-02         -0.318544057474234E-02
  |    2          0.134613400219032E-01          0.285024962984394E-07          0.267807532590122E-07
  |    3          0.459853873574856E-03          0.318543632644722E-02          0.318542786251598E-02
  |    4          0.272039974052440E-03          0.277506938350453E-03          0.277484016527121E-03
  |    5          0.171391834685395E-03         -0.402370692997409E-08          0.128918924994038E-08
  |    6          0.270520422955469E-03         -0.277664475449086E-03         -0.277673180482281E-03
  |    7          0.144473023206898E-01          0.512713272446789E-04          0.109449022059541E-07
  |    8          0.144482814080287E-01          0.396650584237802E-06         -0.509021313935769E-04
  |    9          0.109814915677535E-02         -0.267654083798097E-05          0.420731282893621E-06
  |   10          0.466918963126392E-03         -0.876295533278206E-06         -0.330924619467759E-06
  |   11          0.466512072866561E-03         -0.327286749879462E-06          0.806961213765010E-06
  |   12          0.109712304541722E-02         -0.493306032477536E-06          0.299063760619170E-05
  |   13          0.456253095224984E-03          0.321502709317041E-02         -0.321462818631092E-02
  |   14          0.894934004183427E-03          0.164164098163799E-06         -0.291543246121340E-02
  |   15          0.290606915238854E-03         -0.100230285649372E-03         -0.159541478197569E-06
  |   16          0.272625131015789E-03         -0.280474078591251E-03          0.280296164388911E-03
  |   17          0.291137984288641E-03         -0.153480309357925E-06          0.100422077604930E-03
  |   18          0.896452706744581E-03          0.291453918884263E-02          0.634859729245536E-06
  |   19         -0.588264441716524E-03          0.459789378484797E-03         -0.643808870798734E-03
  |   20         -0.301209237476018E-03          0.136949316630798E-07         -0.574181426536628E-07
  |   21         -0.588703015599948E-03         -0.459768063289049E-03          0.643846072038213E-03
  |   22         -0.588287806861018E-03         -0.643814391413281E-03          0.459793339380723E-03
  |   23         -0.301201681523917E-03         -0.631524165717553E-07          0.201761033536507E-07
  |   24         -0.588691185091410E-03          0.643836436012666E-03         -0.459772761203705E-03
  |   25         -0.575290050552067E-03          0.106828268398685E-05          0.917982695278538E-08
  |   26         -0.127252724647063E-02         -0.532740231512205E-06          0.243544763641080E-05
  |   27         -0.170561195753310E-02         -0.190307459063154E-05          0.185719094014971E-05
  |   28         -0.127405316608035E-02         -0.187016780499191E-05         -0.212325457661516E-06
  |   29         -0.575797788510007E-03          0.145199178973404E-06         -0.120444463333785E-05
  |   30         -0.520318258219003E-03          0.133553664440699E-05         -0.114500449636850E-05
  |   31         -0.587457016753397E-03         -0.455865886523464E-03         -0.638569270127254E-03
  |   32         -0.339633472381555E-02          0.208735377130889E-06         -0.734543566527732E-03
  |   33         -0.339560090181064E-02          0.734294124130646E-03         -0.394704722552815E-06
  |   34         -0.586474238993573E-03          0.638350997827413E-03          0.455504979940412E-03
  |   35         -0.222219628126822E-03          0.298406249379692E-06          0.119518802733161E-03
  |   36         -0.222618454432870E-03         -0.119136998124818E-03          0.291208298125851E-06
  |   37          0.144473042342790E-01          0.161326657901683E-07          0.512694521040663E-04
  |   38          0.144482695870016E-01         -0.508977173254062E-04          0.394868298331683E-06
  |   39          0.109813222571922E-02          0.444265172540015E-06         -0.267614403872161E-05
  |   40          0.466910517452128E-03         -0.339951274227179E-06         -0.872039970152964E-06
  |   41          0.466502936429949E-03          0.795743269445126E-06         -0.318600460263816E-06
  |   42          0.109710835412258E-02          0.298756230482044E-05         -0.490276046636735E-06
  |   43          0.154294424145569E-01          0.888157679185431E-06          0.893699520476667E-06
  |   44          0.455726366849398E-03         -0.321489065379855E-02         -0.321487638613739E-02
  |   45          0.272566280204026E-03         -0.280487101342755E-03         -0.280475542953627E-03
  |   46          0.171889911341001E-03         -0.191061505849799E-06         -0.156069223333008E-06
  |   47          0.274143246448820E-03          0.279861917938715E-03          0.279868217782876E-03
  |   48          0.457268789978671E-03          0.321477000985197E-02          0.321477969039683E-02
  |   49         -0.341239653068126E-02          0.381759795045390E-08         -0.744498310291970E-03
  |   50         -0.589020779672441E-03         -0.459799692907623E-03         -0.643872093760844E-03
  |   51         -0.223340927179136E-03         -0.121324782168154E-03         -0.144940334585945E-07
  |   52         -0.222983213539760E-03         -0.186556318567450E-07          0.121341043405100E-03
  |   53         -0.589200924746521E-03          0.643762882189527E-03          0.459800966602707E-03
  |   54         -0.341263644784457E-02          0.744552333131172E-03          0.664584051122576E-07
  |   55         -0.127404441008613E-02         -0.239627244500651E-06         -0.188783952084066E-05
  |   56         -0.575817137115136E-03         -0.118956540857015E-05          0.169399710056797E-06
  |   57         -0.520340214904260E-03         -0.112114592677995E-05          0.133903203937005E-05
  |   58         -0.575292948141984E-03          0.864902458728192E-08          0.106673185020229E-05
  |   59         -0.127251861163288E-02          0.243177824847227E-05         -0.527738692850098E-06
  |   60         -0.170559743149509E-02          0.184843028478457E-05         -0.190024294797624E-05
  |   61         -0.300746156594190E-03         -0.377409616710583E-06          0.233517141966712E-06
  |   62         -0.586859768126278E-03          0.455587578288686E-03         -0.637963666883416E-03
  |   63         -0.586796249934771E-03          0.638364967630213E-03         -0.455838286894648E-03
  |   64         -0.300734539469824E-03          0.227096711938826E-06         -0.367038729185743E-06
  |   65         -0.586854243439158E-03         -0.637971517413950E-03          0.455587921493768E-03
  |   66         -0.586792899291589E-03         -0.455836881962744E-03          0.638369685970110E-03
  |   67          0.892216729971144E-03          0.216447396149238E-07         -0.289863709002213E-02
  |   68          0.459453348166537E-03          0.318546496898376E-02         -0.318543689788067E-02
  |   69          0.893367442719939E-03          0.289860742967257E-02          0.571276457616680E-08
  |   70          0.290541325626899E-03         -0.234474922440245E-07          0.982609499897220E-04
  |   71          0.271680333993585E-03         -0.277665509393422E-03          0.277510828114160E-03
  |   72          0.290149406223277E-03         -0.982454307247539E-04         -0.192462267327506E-07
  |   73          0.456245653749843E-03         -0.321462044967275E-02          0.321503050432887E-02
  |   74          0.894930725297540E-03         -0.291543677661177E-02          0.148617867743933E-06
  |   75          0.290632087659473E-03         -0.162293063122877E-06         -0.100256721769128E-03
  |   76          0.272650140178998E-03          0.280314515013461E-03         -0.280474923356897E-03
  |   77          0.291125781668046E-03          0.100426987897682E-03         -0.145987111937853E-06
  |   78          0.896435993228713E-03          0.638869122501257E-06          0.291454585149614E-02
  |   79         -0.341240884645496E-02         -0.744496196764947E-03          0.415080021564938E-08
  |   80         -0.588974953295497E-03         -0.643845932136650E-03         -0.459802462788622E-03
  |   81         -0.223295788518652E-03         -0.173976921487962E-07         -0.121322167614127E-03
  |   82         -0.222988063169263E-03          0.121334931944025E-03         -0.109493616582613E-07
  |   83         -0.589192093497825E-03          0.459795338268827E-03          0.643767201114741E-03
  |   84         -0.341263031335884E-02          0.613372587892055E-07          0.744558660355121E-03
  |   85         -0.170614565284495E-02         -0.201448043031923E-05         -0.204198214223519E-05
  |   86         -0.520401315490323E-03         -0.111736700397538E-05         -0.115560037331532E-05
  |   87         -0.324132074810265E-03          0.154047471252221E-06          0.160607464480014E-06
  |   88         -0.518794822788821E-03          0.168679425778024E-05          0.168130998785210E-05
  |   89         -0.170501800839224E-02          0.189179595380746E-05          0.188834038567376E-05
  |   90         -0.287512730917409E-02         -0.557895395067529E-06         -0.569085478106272E-06
  |   91         -0.586453752694832E-03          0.455521131420123E-03          0.638323244853592E-03
  |   92         -0.222234773158444E-03          0.119539214129715E-03          0.302170528201675E-06
  |   93         -0.222608101411307E-03          0.290015810148143E-06         -0.119139538175142E-03
  |   94         -0.587464846937578E-03         -0.638570854791936E-03         -0.455863439134092E-03
  |   95         -0.339633260609773E-02         -0.734544675975996E-03          0.214174128620412E-06
  |   96         -0.339559029956301E-02         -0.391476555948100E-06          0.734288746009631E-03
  |   97          0.290535745303439E-03          0.982663370210020E-04         -0.129179154768023E-07
  |   98          0.271676520447594E-03          0.277524147491618E-03         -0.277662740726568E-03
  |   99          0.290163891872399E-03         -0.372397048485390E-07         -0.982458059306634E-04
  |  100          0.892224010046337E-03         -0.289864262553195E-02          0.337806946448104E-07
  |  101          0.459453179490755E-03         -0.318543656843699E-02          0.318546377480913E-02
  |  102          0.893373497868674E-03         -0.101963288570795E-08          0.289860561819428E-02
  |  103          0.109778643664007E-02         -0.265295068264552E-05         -0.840513781563492E-06
  |  104          0.109786832349238E-02          0.489444473722388E-06          0.298132325822035E-05
  |  105          0.669586708515371E-03         -0.378048188652555E-06          0.397672220773975E-06
  |  106          0.109779299499771E-02         -0.841291739966698E-06         -0.265716198568860E-05
  |  107          0.109787538267595E-02          0.298380869318403E-05          0.487820894564753E-06
  |  108          0.669576728426890E-03          0.400017502412107E-06         -0.370410103027910E-06
  |  109         -0.117622498311294E+00         -0.664367140361455E-06         -0.663005426891026E-06
  |  110         -0.477497155947050E-02         -0.972071387263576E-05         -0.971004899483832E-05
  |  111         -0.687155592053156E-03          0.203683456636944E-06          0.215869410770613E-06
  |  112         -0.398426817326211E-03         -0.275134726871634E-06         -0.268185634185594E-06
  |  113         -0.688662882108337E-03          0.510400059052583E-06          0.507333654826336E-06
  |  114         -0.477629363233544E-02          0.994024286354681E-05          0.993525149022965E-05
  |  115          0.625994626217320E-02          0.201977770360077E-06          0.105329707171903E-01
  |  116         -0.564368707643672E-03          0.863840642574331E-03          0.923004455663203E-03
  |  117         -0.216660899261753E-03          0.220217647346112E-03         -0.800409950301444E-07
  |  118         -0.216389262895209E-03         -0.304907931455811E-06         -0.220351670972197E-03
  |  119         -0.563676679378602E-03         -0.924046489053123E-03         -0.863839065792275E-03
  |  120          0.626026448623760E-02         -0.105326921653655E-01          0.108289857141474E-05
  |  121          0.144139796503158E-02         -0.712529956665577E-08          0.127959769205317E-02
  |  122          0.380944459571124E-03          0.220902709027378E-03         -0.102940195406218E-07
  |  123          0.358496855084066E-03          0.243164235025138E-03         -0.243166961800265E-03
  |  124          0.381383993779777E-03         -0.522115986539943E-07         -0.220887172600157E-03
  |  125          0.144160003267179E-02         -0.127968848309736E-02          0.396829196026443E-07
  |  126          0.179776977797564E-02         -0.927637688282142E-03          0.927746336078750E-03
  |  127          0.101709535070113E-02         -0.120483775223420E-06          0.209833955769639E-06
  |  128          0.101115800402985E-02          0.124679130528018E-05         -0.188135849910508E-05
  |  129          0.101089995778123E-02          0.199788054789377E-05         -0.145991087371772E-05
  |  130          0.101709747416117E-02          0.202202185059734E-06         -0.102482020954159E-06
  |  131          0.101116074271863E-02         -0.188221371498272E-05          0.124914797514246E-05
  |  132          0.101090652400258E-02         -0.146275887734915E-05          0.200093432604358E-05
  |  133          0.143364147377876E-02          0.851454089747721E-06          0.126901724364914E-02
  |  134          0.178852706426026E-02         -0.919062376850606E-03          0.919768303541647E-03
  |  135          0.143509085810235E-02         -0.126990330695482E-02          0.121948501781215E-07
  |  136          0.380726153956503E-03          0.225087311371931E-06         -0.218489804637813E-03
  |  137          0.358039813984558E-03          0.240207108022848E-03         -0.240213755568931E-03
  |  138          0.379972774400199E-03          0.218252527752959E-03         -0.681571986715335E-06
  |  139          0.620936721945111E-02          0.115753490015893E-06          0.104276299089732E-01
  |  140          0.620771213718679E-02         -0.104275827009383E-01          0.287273381181586E-07
  |  141         -0.562700232686120E-03         -0.915410239454151E-03         -0.855917664504348E-03
  |  142         -0.216369548280289E-03         -0.103463216142643E-07         -0.217474192061448E-03
  |  143         -0.216139547982291E-03          0.217466965826543E-03         -0.238461496335864E-07
  |  144         -0.560682617433280E-03          0.855906841117218E-03          0.915271606427075E-03
  |  145          0.625994334850098E-02          0.105329627255592E-01          0.195028006257465E-06
  |  146         -0.564369571276875E-03          0.922991957056100E-03          0.863841966241134E-03
  |  147         -0.216656005820041E-03         -0.827833858766862E-07          0.220231024617423E-03
  |  148         -0.216388013107606E-03         -0.220344664242793E-03         -0.296062672890433E-06
  |  149         -0.563676755515653E-03         -0.863839462735056E-03         -0.924044906396534E-03
  |  150          0.626026614378445E-02          0.108409875241132E-05         -0.105326914260221E-01
  |  151          0.179728891876416E-02          0.927648631879311E-03          0.927646009264489E-03
  |  152          0.358325038383230E-03          0.243163637159684E-03          0.243167367167598E-03
  |  153          0.208822398502267E-03         -0.171786741173149E-07         -0.865273086376036E-08
  |  154          0.358633695668161E-03         -0.243197977999214E-03         -0.243196319532943E-03
  |  155          0.179796757591748E-02         -0.927692267963255E-03         -0.927690158484803E-03
  |  156          0.359245890597778E-02          0.904170159377061E-07          0.905563962787742E-07
  |  157          0.101090219086094E-02          0.140013530943935E-05          0.197898162608690E-05
  |  158          0.454881259797747E-03          0.940544571263206E-06          0.695295299831384E-07
  |  159          0.455151534545645E-03          0.227507346772440E-06         -0.823390014001975E-06
  |  160          0.101146562882286E-02         -0.118136360605041E-05         -0.147149405252212E-05
  |  161          0.194405425208875E-02         -0.320359953689863E-05         -0.817355620624761E-06
  |  162          0.194376557940740E-02         -0.110011291100146E-06          0.299223205216356E-05
  |  163          0.380728959688278E-03         -0.218487388020064E-03          0.235530051801582E-06
  |  164          0.358040593579350E-03         -0.240212468101363E-03          0.240207055644850E-03
  |  165          0.379977034662208E-03         -0.683590230632448E-06          0.218252917639235E-03
  |  166          0.143363553186894E-02          0.126901776493554E-02          0.859701474741339E-06
  |  167          0.178852912183490E-02          0.919766371041043E-03         -0.919059111326169E-03
  |  168          0.143510179362493E-02         -0.636788433362690E-08         -0.126990241163678E-02
  |  169         -0.562462211984102E-03         -0.855876389196044E-03          0.915378579925178E-03
  |  170         -0.562602112860941E-03         -0.915392754394149E-03          0.855922508606063E-03
  |  171         -0.461636820705358E-03         -0.540306929539059E-07          0.172379943925665E-07
  |  172         -0.562466567543203E-03          0.915378158405711E-03         -0.855873767104802E-03
  |  173         -0.562604745705413E-03          0.855923817852898E-03         -0.915389701960565E-03
  |  174         -0.461635298109988E-03          0.123381223753715E-07         -0.481965255484627E-07
  |  175         -0.477572740729509E-02         -0.105286887706655E-04          0.998509161340432E-05
  |  176         -0.215556559602345E-02          0.426991531356180E-05          0.559008458834997E-09
  |  177         -0.761985046702859E-03         -0.274628896015975E-06          0.163113852227147E-05
  |  178         -0.687233090922129E-03         -0.214702606196903E-06          0.217664160986218E-06
  |  179         -0.762738496436699E-03         -0.135933108564274E-05          0.593399057508970E-06
  |  180         -0.215706504611191E-02         -0.720971632707438E-06         -0.359889011044017E-05
  |  181          0.144140914845699E-02          0.127959555089612E-02         -0.209219143354277E-08
  |  182          0.380971908848995E-03         -0.134382591373036E-07          0.220901086268380E-03
  |  183          0.358509932725897E-03         -0.243161060026819E-03          0.243163220159949E-03
  |  184          0.381385133078248E-03         -0.220881133214432E-03         -0.509529789130016E-07
  |  185          0.144159847502907E-02          0.372809819014852E-07         -0.127968817605308E-02
  |  186          0.179776605230957E-02          0.927744011641238E-03         -0.927636582058519E-03
  |  187          0.101089838820282E-02          0.198920554229172E-05          0.138783468635648E-05
  |  188          0.454890606270504E-03          0.675709229034875E-07          0.931570698890392E-06
  |  189          0.455150316188485E-03         -0.821065938691415E-06          0.225940503563809E-06
  |  190          0.101146665742501E-02         -0.147156334736095E-05         -0.118329729730058E-05
  |  191          0.194405783185670E-02         -0.815733443750411E-06         -0.320146668229700E-05
  |  192          0.194376060801805E-02          0.300137075266659E-05         -0.105349726838950E-06
  |  193          0.358145188494597E-03         -0.240212183124391E-03         -0.240224866264834E-03
  |  194          0.208943726278731E-03          0.243954814691417E-06          0.241743510764434E-06
  |  195          0.356614037748766E-03          0.239359876241525E-03          0.239355984294596E-03
  |  196          0.178739226323200E-02          0.919760135420444E-03          0.919756316197189E-03
  |  197          0.356649586407693E-02          0.861982450678132E-06          0.861057659643851E-06
  |  198          0.178930153220864E-02         -0.919987878568532E-03         -0.920000406005163E-03
  |  199         -0.216370872941733E-03         -0.217472567890395E-03         -0.949424326199526E-08
  |  200         -0.216137616697596E-03         -0.295191673998805E-07          0.217464001050694E-03
  |  201         -0.560682815935795E-03          0.915276054673733E-03          0.855901528273864E-03
  |  202          0.620936529455318E-02          0.104276336259101E-01          0.120087090742335E-06
  |  203          0.620771251831169E-02          0.288496869455245E-07         -0.104275843935648E-01
  |  204         -0.562692357941687E-03         -0.855918943116583E-03         -0.915410068103591E-03
  |  205         -0.687228541997606E-03          0.216011189107877E-06         -0.218878399862070E-06
  |  206         -0.762724295680055E-03          0.586343023430915E-06         -0.136197632887891E-05
  |  207         -0.215706319723358E-02         -0.359896723279682E-05         -0.717788284305660E-06
  |  208         -0.477572590280627E-02          0.998763268510391E-05         -0.105281039360233E-04
  |  209         -0.215555532834078E-02          0.133499949655158E-07          0.426966618343572E-05
  |  210         -0.761982986560107E-03          0.163386637915067E-05         -0.278871031130162E-06
  |  211         -0.564366583598966E-03          0.922987702499419E-03         -0.863767864339182E-03
  |  212         -0.461299507026297E-03         -0.294703024374263E-06          0.133068453238418E-06
  |  213         -0.564090780797397E-03         -0.923189542958452E-03          0.864126778108754E-03
  |  214         -0.564361418952570E-03         -0.863767807118678E-03          0.922989596225672E-03
  |  215         -0.461303661438938E-03          0.139754494580330E-06         -0.297819957127254E-06
  |  216         -0.564096757944293E-03          0.864132586867956E-03         -0.923195060638363E-03

  Removing unitary transformations (pure translations, rotations) from forces on atoms.
  Atomic forces before filtering:
  | Net force on center of mass :  -0.620240E-16  0.172482E-16 -0.682961E-16 eV/A
  Atomic forces after filtering:
  | Net force on center of mass :  -0.135024E-15  0.170740E-16 -0.681219E-16 eV/A
 ------------------------------------------------------------------------
 Atomic structure that was used in the preceding time step of the wrapper
                         x [A]             y [A]             z [A]
  lattice_vector        12.67881443        0.00000000        0.00000000
  lattice_vector         0.00000000       12.67881443        0.00000000
  lattice_vector        -0.00000000       -0.00000000       12.67881443

            atom        -0.00000000        0.00000000       -0.00000000  Mg
            atom        -0.00000000        2.11313574        2.11313574  Mg
            atom        -0.00000000        4.22627148        4.22627148  Mg
            atom        -0.00000000        6.33940721        6.33940721  Mg
            atom        -0.00000000       -4.22627148       -4.22627148  Mg
            atom        -0.00000000       -2.11313574       -2.11313574  Mg
            atom         2.11313574        0.00000000        2.11313574  Mg
            atom         2.11313574        2.11313574        4.22627148  Mg
            atom         2.11313574        4.22627148        6.33940721  Mg
            atom         2.11313574        6.33940721       -4.22627148  Mg
            atom         2.11313574       -4.22627148       -2.11313574  Mg
            atom         2.11313574       -2.11313574       -0.00000000  Mg
            atom         4.22627148        0.00000000        4.22627148  Mg
            atom         4.22627148        2.11313574        6.33940721  Mg
            atom         4.22627148        4.22627148       -4.22627148  Mg
            atom         4.22627148        6.33940721       -2.11313574  Mg
            atom         4.22627148       -4.22627148       -0.00000000  Mg
            atom         4.22627148       -2.11313574        2.11313574  Mg
            atom         6.33940721        0.00000000        6.33940721  Mg
            atom         6.33940721        2.11313574       -4.22627148  Mg
            atom         6.33940721        4.22627148       -2.11313574  Mg
            atom         6.33940721        6.33940721        0.00000000  Mg
            atom         6.33940721       -4.22627148        2.11313574  Mg
            atom         6.33940721       -2.11313574        4.22627148  Mg
            atom        -4.22627148       -0.00000000       -4.22627148  Mg
            atom        -4.22627148        2.11313574       -2.11313574  Mg
            atom        -4.22627148        4.22627148        0.00000000  Mg
            atom        -4.22627148        6.33940721        2.11313574  Mg
            atom        -4.22627148       -4.22627148        4.22627148  Mg
            atom        -4.22627148       -2.11313574        6.33940721  Mg
            atom        -2.11313574       -0.00000000       -2.11313574  Mg
            atom        -2.11313574        2.11313574        0.00000000  Mg
            atom        -2.11313574        4.22627148        2.11313574  Mg
            atom        -2.11313574        6.33940721        4.22627148  Mg
            atom        -2.11313574       -4.22627148        6.33940721  Mg
            atom        -2.11313574       -2.11313574       -4.22627148  Mg
            atom         2.11313574        2.11313574        0.00000000  Mg
            atom         2.11313574        4.22627148        2.11313574  Mg
            atom         2.11313574        6.33940721        4.22627148  Mg
            atom         2.11313574       -4.22627148        6.33940721  Mg
            atom         2.11313574       -2.11313574       -4.22627148  Mg
            atom         2.11313574       -0.00000000       -2.11313574  Mg
            atom         4.22627148        2.11313574        2.11313574  Mg
            atom         4.22627148        4.22627148        4.22627148  Mg
            atom         4.22627148        6.33940721        6.33940721  Mg
            atom         4.22627148       -4.22627148       -4.22627148  Mg
            atom         4.22627148       -2.11313574       -2.11313574  Mg
            atom         4.22627148        0.00000000       -0.00000000  Mg
            atom         6.33940721        2.11313574        4.22627148  Mg
            atom         6.33940721        4.22627148        6.33940721  Mg
            atom         6.33940721        6.33940721       -4.22627148  Mg
            atom         6.33940721       -4.22627148       -2.11313574  Mg
            atom         6.33940721       -2.11313574       -0.00000000  Mg
            atom         6.33940721        0.00000000        2.11313574  Mg
            atom        -4.22627148        2.11313574        6.33940721  Mg
            atom        -4.22627148        4.22627148       -4.22627148  Mg
            atom        -4.22627148        6.33940721       -2.11313574  Mg
            atom        -4.22627148       -4.22627148       -0.00000000  Mg
            atom        -4.22627148       -2.11313574        2.11313574  Mg
            atom        -4.22627148       -0.00000000        4.22627148  Mg
            atom        -2.11313574        2.11313574       -4.22627148  Mg
            atom        -2.11313574        4.22627148       -2.11313574  Mg
            atom        -2.11313574        6.33940721        0.00000000  Mg
            atom        -2.11313574       -4.22627148        2.11313574  Mg
            atom        -2.11313574       -2.11313574        4.22627148  Mg
            atom        -2.11313574       -0.00000000        6.33940721  Mg
            atom         0.00000000        2.11313574       -2.11313574  Mg
            atom         0.00000000        4.22627148        0.00000000  Mg
            atom         0.00000000        6.33940721        2.11313574  Mg
            atom        -0.00000000       -4.22627148        4.22627148  Mg
            atom        -0.00000000       -2.11313574        6.33940721  Mg
            atom         0.00000000       -0.00000000       -4.22627148  Mg
            atom         4.22627148        4.22627148        0.00000000  Mg
            atom         4.22627148        6.33940721        2.11313574  Mg
            atom         4.22627148       -4.22627148        4.22627148  Mg
            atom         4.22627148       -2.11313574        6.33940721  Mg
            atom         4.22627148       -0.00000000       -4.22627148  Mg
            atom         4.22627148        2.11313574       -2.11313574  Mg
            atom         6.33940721        4.22627148        2.11313574  Mg
            atom         6.33940721        6.33940721        4.22627148  Mg
            atom         6.33940721       -4.22627148        6.33940721  Mg
            atom         6.33940721       -2.11313574       -4.22627148  Mg
            atom         6.33940721       -0.00000000       -2.11313574  Mg
            atom         6.33940721        2.11313574        0.00000000  Mg
            atom        -4.22627148        4.22627148        4.22627148  Mg
            atom        -4.22627148        6.33940721        6.33940721  Mg
            atom        -4.22627148       -4.22627148       -4.22627148  Mg
            atom        -4.22627148       -2.11313574       -2.11313574  Mg
            atom        -4.22627148       -0.00000000        0.00000000  Mg
            atom        -4.22627148        2.11313574        2.11313574  Mg
            atom        -2.11313574        4.22627148        6.33940721  Mg
            atom        -2.11313574        6.33940721       -4.22627148  Mg
            atom        -2.11313574       -4.22627148       -2.11313574  Mg
            atom        -2.11313574       -2.11313574       -0.00000000  Mg
            atom        -2.11313574       -0.00000000        2.11313574  Mg
            atom        -2.11313574        2.11313574        4.22627148  Mg
            atom         0.00000000        4.22627148       -4.22627148  Mg
            atom         0.00000000        6.33940721       -2.11313574  Mg
            atom        -0.00000000       -4.22627148       -0.00000000  Mg
            atom        -0.00000000       -2.11313574        2.11313574  Mg
            atom        -0.00000000       -0.00000000        4.22627148  Mg
            atom        -0.00000000        2.11313574        6.33940721  Mg
            atom         2.11313574        4.22627148       -2.11313574  Mg
            atom         2.11313574        6.33940721        0.00000000  Mg
            atom         2.11313574       -4.22627148        2.11313574  Mg
            atom         2.11313574       -2.11313574        4.22627148  Mg
            atom         2.11313574       -0.00000000        6.33940721  Mg
            atom         2.11313574        2.11313574       -4.22627148  Mg
            atom         2.12313574        2.11313574        2.11313574  O
            atom         2.11313574        4.22627148        4.22627148  O
            atom         2.11313574        6.33940721        6.33940721  O
            atom         2.11313574       -4.22627148       -4.22627148  O
            atom         2.11313574       -2.11313574       -2.11313574  O
            atom         2.11313574        0.00000000       -0.00000000  O
            atom         4.22627148        2.11313574        4.22627148  O
            atom         4.22627148        4.22627148        6.33940721  O
            atom         4.22627148        6.33940721       -4.22627148  O
            atom         4.22627148       -4.22627148       -2.11313574  O
            atom         4.22627148       -2.11313574       -0.00000000  O
            atom         4.22627148        0.00000000        2.11313574  O
            atom         6.33940721        2.11313574        6.33940721  O
            atom         6.33940721        4.22627148       -4.22627148  O
            atom         6.33940721        6.33940721       -2.11313574  O
            atom         6.33940721       -4.22627148       -0.00000000  O
            atom         6.33940721       -2.11313574        2.11313574  O
            atom         6.33940721        0.00000000        4.22627148  O
            atom        -4.22627148        2.11313574       -4.22627148  O
            atom        -4.22627148        4.22627148       -2.11313574  O
            atom        -4.22627148        6.33940721        0.00000000  O
            atom        -4.22627148       -4.22627148        2.11313574  O
            atom        -4.22627148       -2.11313574        4.22627148  O
            atom        -4.22627148        0.00000000        6.33940721  O
            atom        -2.11313574        2.11313574       -2.11313574  O
            atom        -2.11313574        4.22627148        0.00000000  O
            atom        -2.11313574        6.33940721        2.11313574  O
            atom        -2.11313574       -4.22627148        4.22627148  O
            atom        -2.11313574       -2.11313574        6.33940721  O
            atom        -2.11313574        0.00000000       -4.22627148  O
            atom        -0.00000000        2.11313574       -0.00000000  O
            atom        -0.00000000        4.22627148        2.11313574  O
            atom        -0.00000000        6.33940721        4.22627148  O
            atom        -0.00000000       -4.22627148        6.33940721  O
            atom        -0.00000000       -2.11313574       -4.22627148  O
            atom        -0.00000000        0.00000000       -2.11313574  O
            atom         4.22627148        4.22627148        2.11313574  O
            atom         4.22627148        6.33940721        4.22627148  O
            atom         4.22627148       -4.22627148        6.33940721  O
            atom         4.22627148       -2.11313574       -4.22627148  O
            atom         4.22627148       -0.00000000       -2.11313574  O
            atom         4.22627148        2.11313574       -0.00000000  O
            atom         6.33940721        4.22627148        4.22627148  O
            atom         6.33940721        6.33940721        6.33940721  O
            atom         6.33940721       -4.22627148       -4.22627148  O
            atom         6.33940721       -2.11313574       -2.11313574  O
            atom         6.33940721        0.00000000       -0.00000000  O
            atom         6.33940721        2.11313574        2.11313574  O
            atom        -4.22627148        4.22627148        6.33940721  O
            atom        -4.22627148        6.33940721       -4.22627148  O
            atom        -4.22627148       -4.22627148       -2.11313574  O
            atom        -4.22627148       -2.11313574       -0.00000000  O
            atom        -4.22627148       -0.00000000        2.11313574  O
            atom        -4.22627148        2.11313574        4.22627148  O
            atom        -2.11313574        4.22627148       -4.22627148  O
            atom        -2.11313574        6.33940721       -2.11313574  O
            atom        -2.11313574       -4.22627148       -0.00000000  O
            atom        -2.11313574       -2.11313574        2.11313574  O
            atom        -2.11313574       -0.00000000        4.22627148  O
            atom        -2.11313574        2.11313574        6.33940721  O
            atom         0.00000000        4.22627148       -2.11313574  O
            atom         0.00000000        6.33940721        0.00000000  O
            atom        -0.00000000       -4.22627148        2.11313574  O
            atom        -0.00000000       -2.11313574        4.22627148  O
            atom        -0.00000000       -0.00000000        6.33940721  O
            atom         0.00000000        2.11313574       -4.22627148  O
            atom         2.11313574        4.22627148       -0.00000000  O
            atom         2.11313574        6.33940721        2.11313574  O
            atom         2.11313574       -4.22627148        4.22627148  O
            atom         2.11313574       -2.11313574        6.33940721  O
            atom         2.11313574       -0.00000000       -4.22627148  O
            atom         2.11313574        2.11313574       -2.11313574  O
            atom         6.33940721        6.33940721        2.11313574  O
            atom         6.33940721       -4.22627148        4.22627148  O
            atom         6.33940721       -2.11313574        6.33940721  O
            atom         6.33940721       -0.00000000       -4.22627148  O
            atom         6.33940721        2.11313574       -2.11313574  O
            atom         6.33940721        4.22627148       -0.00000000  O
            atom        -4.22627148        6.33940721        4.22627148  O
            atom        -4.22627148       -4.22627148        6.33940721  O
            atom        -4.22627148       -2.11313574       -4.22627148  O
            atom        -4.22627148       -0.00000000       -2.11313574  O
            atom        -4.22627148        2.11313574        0.00000000  O
            atom        -4.22627148        4.22627148        2.11313574  O
            atom        -2.11313574        6.33940721        6.33940721  O
            atom        -2.11313574       -4.22627148       -4.22627148  O
            atom        -2.11313574       -2.11313574       -2.11313574  O
            atom        -2.11313574       -0.00000000        0.00000000  O
            atom        -2.11313574        2.11313574        2.11313574  O
            atom        -2.11313574        4.22627148        4.22627148  O
            atom         0.00000000        6.33940721       -4.22627148  O
            atom        -0.00000000       -4.22627148       -2.11313574  O
            atom        -0.00000000       -2.11313574       -0.00000000  O
            atom        -0.00000000       -0.00000000        2.11313574  O
            atom        -0.00000000        2.11313574        4.22627148  O
            atom        -0.00000000        4.22627148        6.33940721  O
            atom         2.11313574        6.33940721       -2.11313574  O
            atom         2.11313574       -4.22627148       -0.00000000  O
            atom         2.11313574       -2.11313574        2.11313574  O
            atom         2.11313574       -0.00000000        4.22627148  O
            atom         2.11313574        2.11313574        6.33940721  O
            atom         2.11313574        4.22627148       -4.22627148  O
            atom         4.22627148        6.33940721       -0.00000000  O
            atom         4.22627148       -4.22627148        2.11313574  O
            atom         4.22627148       -2.11313574        4.22627148  O
            atom         4.22627148       -0.00000000        6.33940721  O
            atom         4.22627148        2.11313574       -4.22627148  O
            atom         4.22627148        4.22627148       -2.11313574  O
 ------------------------------------------------------------------------
  @ DRIVER MODE: Message from server: STATUS
  @ DRIVER MODE: Message from server: GETFORCE
  @ DRIVER MODE: Returning v,forces,stress

------------------------------------------------------------------------------
  Final output of selected total energy values:

  The following output summarizes some interesting total energy values
  at the end of a run (AFTER all relaxation, molecular dynamics, etc.).

  | Total energy of the DFT / Hartree-Fock s.c.f. calculation      :        -808886.468076903 eV
  | Final zero-broadening corrected energy (caution - metals only) :        -808886.468076903 eV
  | For reference only, the value of 1 Hartree used in FHI-aims is :             27.211384500 eV
  | For reference only, the overall average (free atom contribution
  | + realspace contribution) of the electrostatic potential after
  | s.c.f. convergence is                                          :            -16.776204048 eV

  Before relying on these values, please be sure to understand exactly which
  total energy value is referred to by a given number. Different objects may
  all carry the same name 'total energy'. Definitions:

  Total energy of the DFT / Hartree-Fock s.c.f. calculation:
  | Note that this energy does not include ANY quantities calculated after the
  | s.c.f. cycle, in particular not ANY RPA, MP2, etc. many-body perturbation terms.

  Final zero-broadening corrected energy:
  | For metallic systems only, a broadening of the occupation numbers at the Fermi
  | level can be extrapolated back to zero broadening by an electron-gas inspired
  | formula. For all systems that are not real metals, this value can be
  | meaningless and should be avoided.

------------------------------------------------------------------------------
  Methods described in the following list of references were used in this FHI-aims run.
  If you publish the results, please make sure to cite these reference if they apply.
  FHI-aims is an academic code, and for our developers (often, Ph.D. students
  and postdocs), scientific credit in the community is essential.
  Thank you for helping us!

  For any use of FHI-aims, please cite:

    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    Xinguo Ren, Karsten Reuter, and Matthias Scheffler
    'Ab initio molecular simulations with numeric atom-centered orbitals'
    Computer Physics Communications 180, 2175-2196 (2009)
    http://dx.doi.org/10.1016/j.cpc.2009.06.022


  The provided symmetry information was generated with SPGlib:

    Atsushi Togo, Yusuke Seto, Dimitar Pashov
    SPGlib 1.7.3 obtained from http://spglib.sourceforge.net
    Copyright (C) 2008 Atsushi Togo


  The ELSI infrastructure was used in your run to solve the Kohn-Sham electronic structure.
  Please check out http://elsi-interchange.org to learn more.
  If scalability is important for your project, please acknowledge ELSI by citing:

    V. W-z. Yu, F. Corsetti, A. Garcia, W. P. Huhn, M. Jacquelin, W. Jia,
    B. Lange, L. Lin, J. Lu, W. Mi, A. Seifitokaldani, A. Vazquez-Mayagoitia,
    C. Yang, H. Yang, and V. Blum
    'ELSI: A unified software interface for Kohn-Sham electronic structure solvers'
    Computer Physics Communications 222, 267-285 (2018).
    http://dx.doi.org/10.1016/j.cpc.2017.09.007


  For the real-space grid partitioning and parallelization used in this calculation, please cite:

    Ville Havu, Volker Blum, Paula Havu, and Matthias Scheffler,
    'Efficient O(N) integration for all-electron electronic structure calculation'
    'using numerically tabulated basis functions'
    Journal of Computational Physics 228, 8367-8379 (2009).
    http://dx.doi.org/10.1016/j.jcp.2009.08.008

  Of course, there are many other important community references, e.g., those cited in the
  above references. Our list is limited to references that describe implementations in the
  FHI-aims code. The reason is purely practical (length of this list) - please credit others as well.

   One or more s.c.f. cycles took more than 50 s.c.f. iterations.
   It might be worth adjusting your convergence settings. For example,
   metals need a larger broadening and a smaller "charge_mix_param" value than semiconductors.

------------------------------------------------------------
          Leaving FHI-aims.
          Date     :  20191125, Time     :  162630.319

          Computational steps:
          | Number of self-consistency cycles          :          102
          | Number of SCF (re)initializations          :            2

          Detailed time accounting                     :  max(cpu_time)    wall_clock(cpu1)
          | Total time                                 :      817.664 s         818.283 s
          | Preparation time                           :        0.217 s           0.415 s
          | Boundary condition initalization           :        3.013 s           3.016 s
          | Grid partitioning                          :        1.198 s           1.178 s
          | Preloading free-atom quantities on grid    :        1.546 s           1.549 s
          | Free-atom superposition energy             :        0.461 s           0.463 s
          | Total time for integrations                :      144.405 s         144.370 s
          | Total time for solution of K.-S. equations :      109.706 s         109.735 s
          | Total time for density & force components  :      262.718 s         262.718 s
          | Total time for mixing & preconditioning    :      133.419 s         133.416 s
          | Total time for Hartree multipole update    :        2.632 s           2.637 s
          | Total time for Hartree multipole sum       :      154.433 s         154.443 s
          | Total time for total energy evaluation     :        0.065 s           0.159 s
          | Total time for scaled ZORA corrections     :        0.000 s           0.000 s

          Partial memory accounting:
          | Residual value for overall tracked memory usage across tasks:     0.000000 MB (should be 0.000000 MB)
          | Peak values for overall tracked memory usage:
          |   Minimum:      232.455 MB (on task  19 after allocating temp_ham_ovlp)
          |   Maximum:      232.789 MB (on task   0 after allocating temp_ham_ovlp)
          |   Average:      232.669 MB
          | Largest tracked array allocation:
          |   Minimum:       56.278 MB (overlap_matrix on task   0)
          |   Maximum:       56.278 MB (overlap_matrix on task   0)
          |   Average:       56.278 MB
          Note:  These values currently only include a subset of arrays which are explicitly tracked.
          The "true" memory usage will be greater.

          Have a nice day.
------------------------------------------------------------
"""

    with open("aims.out", "w") as f:
        f.write(output)
