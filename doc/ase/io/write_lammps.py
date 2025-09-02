from ase_koopmans.io.opls import OPLSff, OPLSStructure

s = OPLSStructure('172_ext.xyz')
opls = OPLSff('172_defs.par')
opls.write_lammps(s, prefix='lmp')
