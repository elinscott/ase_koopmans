def test_ts09():
    from ase import io
    from ase.calculators.vdwcorrection import vdWTkatchenko09prl
    from ase.calculators.emt import EMT
    from ase.build import bulk

    # fake objects for the test
    class FakeHirshfeldPartitioning:
        def __init__(self, calculator):
            self.calculator = calculator

        def initialize(self):
            pass

        def get_effective_volume_ratios(self):
            return [1]

        def get_calculator(self):
            return self.calculator

    class FakeDFTcalculator(EMT):
        def get_xc_functional(self):
            return 'PBE'

    a = 4.05  # Angstrom lattice spacing
    al = bulk('Al', 'fcc', a=a)

    cc = FakeDFTcalculator()
    hp = FakeHirshfeldPartitioning(cc)
    c = vdWTkatchenko09prl(hp, [3])
    al.calc = c
    al.get_potential_energy()

    fname = 'out.traj'
    al.write(fname)

    # check that the output exists
    io.read(fname)
    # maybe assert something about what we just read?

    p = io.read(fname).calc.parameters
    p['calculator']
    p['xc']
    p['uncorrected_energy']
