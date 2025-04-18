from ase_koopmans.build import bulk
from ase_koopmans.calculators.nwchem import NWChem
from ase_koopmans.spectrum.band_structure import calculate_band_structure


def test_bands():
    atoms = bulk('Si')
    path = atoms.cell.bandpath('GXWK', density=10)
    atoms.calc = NWChem(kpts=[2, 2, 2])
    bs = calculate_band_structure(atoms, path)
    print(bs)
    bs.write('bs.json')
