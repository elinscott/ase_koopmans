import pytest
from ase.build import bulk
from ase.dft.bandgap import bandgap
from ase.calculators.calculator import get_calculator_class

kpts = (4, 4, 4)

required = {'abinit': dict(ecut=200, toldfe=0.0001, chksymbreak=0),
            'aims': dict(sc_accuracy_rho=5.e-4,
                         output=['k_point_list',
                                 'k_eigenvalue,' +
                                 str(kpts[0] * kpts[1] * kpts[2])]),
            'elk': dict(tasks=0, rgkmax=5.0),
            'gpaw': dict(mode='pw')}

# XXX This test is in shambles, it's broken for several calculators.
# Not sure whether it works for *any* calculator.
# Someone should rehabilitate this, or we delete it.  --askhl

@pytest.mark.parametrize('name', ['aims', 'elk', 'openmx'])
def test_bandgap(name):
    Calculator = get_calculator_class(name)
    par = required.get(name, {})
    calc = Calculator(label=name + '_bandgap', xc='PBE',
                      # abinit, aims, elk - do not recognize the syntax below:
                      # kpts={'size': kpts, 'gamma': True}, **par)
                      kpts=kpts, **par)
    si = bulk('Si', crystalstructure='diamond', a=5.43)
    si.calc = calc
    si.get_potential_energy()
    print(name, bandgap(si.calc))
    del si.calc
    # test spin-polarization
    calc = Calculator(label=name + '_bandgap_spinpol', xc='PBE',
                      # abinit, aims, elk - do not recognize the syntax below:
                      # kpts={'size': kpts, 'gamma': True}, **par)
                      kpts=kpts, **par)
    si.set_initial_magnetic_moments([-0.1, 0.1])
    # this should not be necessary in the new ase interface standard ...
    if si.get_initial_magnetic_moments().any():  # spin-polarization
        if name == 'aims':
            calc.set(spin='collinear')
        if name == 'elk':
            calc.set(spinpol=True)
    si.set_calculator(calc)
    si.get_potential_energy()
    print(name, bandgap(si.calc))
