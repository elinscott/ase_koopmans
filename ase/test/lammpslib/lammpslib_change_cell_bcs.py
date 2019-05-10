# test that a change in unit cell boundary conditions is
# handled correctly by lammpslib
from __future__ import print_function
import os
from ase.calculators.lammpslib import LAMMPSlib
from ase.lattice.cubic import FaceCenteredCubic

potential_path = os.environ.get('LAMMPS_POTENTIALS_PATH', '.')
cmds = ["pair_style eam/alloy",
        "pair_coeff * * {path}/NiAlH_jea.eam.alloy Ni H"
        "".format(path=potential_path)]
lammps = LAMMPSlib(lmpcmds=cmds,
                   atom_types={'Ni': 1, 'H': 2},
                   log_file='test.log', keep_alive=True)
atoms = FaceCenteredCubic(size=(2,2,2), latticeconstant=3.52, symbol="Ni", pbc=True)
atoms.set_calculator(lammps)

energy_ppp_ref = -142.400000403
energy_ppp = atoms.get_potential_energy()
print("Computed energy with boundary ppp = {}".format(energy_ppp))
diff_ppp = abs((energy_ppp - energy_ppp_ref) / energy_ppp_ref)
assert diff_ppp < 1e-10, "Expected "+str(energy_ppp_ref)

atoms.set_pbc((False,False,True))
energy_ssp_ref = -114.524625705
energy_ssp = atoms.get_potential_energy()
print("Computed energy with boundary ssp = {}".format(energy_ssp))
diff_ssp = abs((energy_ssp - energy_ssp_ref) / energy_ssp_ref)
assert diff_ssp < 1e-10, "Expected "+str(energy_ssp_ref)
