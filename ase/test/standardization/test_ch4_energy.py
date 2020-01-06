import pytest
from ase.build import molecule
from ase.calculators.calculator import get_calculator_class
from ase.units import Ry
from ase.utils import workdir


class CalculatorInputs:
    def __init__(self, name, parameters=None):
        self.name = name
        if parameters is None:
            parameters = {}
        self.parameters = parameters

    def __repr__(self):
        cls = type(self)
        return '{}({}, {})'.format(cls.__name__,
                                   self.name, self.parameters)

    def calc(self):
        cls = get_calculator_class(self.name)
        return cls(**self.parameters)


def inputs(name, **parameters):
    return CalculatorInputs(name, parameters)


def _calculate(spec, name):
    atoms = molecule(name)
    atoms.center(vacuum=3.5)
    with workdir('test-{}'.format(name), mkdir=True):
        atoms.calc = spec.calc()
        return atoms.get_potential_energy()


@pytest.mark.parametrize(
    "spec",
    [
        inputs('abinit', ecut=300, chksymbreak=0, toldfe=1e-4),
        inputs('cp2k'),
        # Siesta gets a result of -3.3 eV, not so accurate then.
        # What should we do about this?
        inputs(
            'espresso',
            ecutwfc=300 / Ry,
            pseudopotentials={'C': 'C.pz-kjpaw.UPF',
                              'H': 'H.pz-kjpaw.UPF'},
        ),
        inputs('gpaw', symmetry='off', mode='pw', txt='gpaw.txt',
               mixer={'beta': 0.6}),
        inputs('octopus', stdout="'stdout.log'", stderr="'stderr.log'"),
        inputs('openmx', energy_cutoff=350),
        pytest.param(inputs('siesta'), marks=pytest.mark.xfail),
    ],
    ids=lambda spec: spec.name)
def test_ch4(tmp_path, spec):
    with workdir(tmp_path, mkdir=True):
        e_ch4 = _calculate(spec, 'CH4')
        e_c2h2 = _calculate(spec, 'C2H2')
        e_h2 = _calculate(spec, 'H2')
        energy = e_ch4 - 0.5 * e_c2h2 - 1.5 * e_h2
        print(energy)
        ref_energy = -2.8
        assert abs(energy - ref_energy) < 0.3
