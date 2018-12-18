from __future__ import print_function
from ase.build import niggli_reduce
from ase.calculators.singlepoint import SinglePointCalculator
from asap3 import EMT
from ase.optimize.precon import PreconLBFGS, PreconFIRE
from ase.ga import set_raw_score


def finalize(atoms, energy=None, forces=None, stress=None):
    atoms.wrap()
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces,
                                 stress=stress)
    atoms.set_calculator(calc)
    raw_score = -atoms.get_potential_energy()
    set_raw_score(atoms, raw_score)


def relax_one(atoms, cellbounds=None):
    calc = EMT()

    # Relax the atoms object sequentially by trying to reduce the unit
    # cell after each couple of relaxation steps
    converged = False
    niter = 0
    while not converged and niter < 100:
        if cellbounds is not None:
            if not cellbounds.is_within_bounds(atoms.get_cell()):
                niggli_reduce(atoms)
            if not cellbounds.is_within_bounds(atoms.get_cell()):
                raise RuntimeError('Niggli reduction did not work; aborting')

        try:
            calc = EMT()
            atoms.set_calculator(calc)
            dyn = PreconLBFGS(atoms, variable_cell=True, maxstep=0.2,
                              a_min=1e-3,
                              use_armijo=True, logfile=None, trajectory=None)
            dyn.run(fmax=5e-2, smax=1e-4, steps=50)
        except RuntimeError:
            dyn = PreconFIRE(atoms, variable_cell=True, dt=0.05, maxmove=0.2,
                             use_armijo=False, logfile=None, trajectory=None)
            dyn.run(fmax=5e-2, smax=1e-4, steps=10)
        converged = dyn.converged()
        niter += 1

    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    s = atoms.get_stress()
    finalize(atoms, energy=e, forces=f, stress=s)
    print('in relax', atoms.calc)
