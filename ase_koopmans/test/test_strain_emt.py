def test_strain_emt():
    """This test checks that the StrainFilter works using the default
    built-in EMT calculator."""

    import numpy as np
    from ase_koopmans.constraints import StrainFilter
    from ase_koopmans.optimize.mdmin import MDMin
    from ase_koopmans.calculators.emt import EMT
    from ase_koopmans.build import bulk

    cu = bulk('Cu', 'fcc', a=3.6)

    class EMTPlus(EMT):
        def get_stress(self, atoms):
            return np.zeros(6)

    cu.calc = EMTPlus()
    f = StrainFilter(cu)
    opt = MDMin(f, dt=0.01)
    opt.run(0.1, steps=2)
