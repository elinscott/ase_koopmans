def test_eam_run():
    import numpy as np

    from ase_koopmans.calculators.eam import EAM

    from ase_koopmans.test.eam_pot import Pt_u3
    from ase_koopmans.build import fcc111
    from io import StringIO

    eam = EAM(potential=StringIO(Pt_u3), form='eam', elements=['Pt'])
    slab = fcc111('Pt', size=(4, 4, 2), vacuum=10.0)
    slab.calc = eam

    assert( abs(-164.277599313 - slab.get_potential_energy()) < 1E-8 )
    assert( abs(6.36379627645 - np.linalg.norm(slab.get_forces()))  < 1E-8 )
