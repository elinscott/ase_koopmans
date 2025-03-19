"""Check that reading and writing masses in .con files is consistent."""

from numpy import asarray
import ase_koopmans.lattice.compounds
import ase_koopmans.data
import ase_koopmans.io

def test_eon_masses():
    # Error tolerance.
    TOL = 1e-8

    data = ase_koopmans.lattice.compounds.B2(['Cs', 'Cl'], latticeconstant=4.123,
                                    size=(3, 3, 3))

    m_Cs = ase_koopmans.data.atomic_masses[ase_koopmans.data.atomic_numbers['Cs']]
    m_Cl = ase_koopmans.data.atomic_masses[ase_koopmans.data.atomic_numbers['Cl']]


    con_file = 'pos.con'
    # Write and read the .con file.
    ase_koopmans.io.write(con_file, data, format='eon')
    data2 = ase_koopmans.io.read(con_file, format='eon')
    # Check masses.
    symbols = asarray(data2.get_chemical_symbols())
    masses = asarray(data2.get_masses())
    assert (abs(masses[symbols == 'Cs'] - m_Cs)).sum() < TOL
    assert (abs(masses[symbols == 'Cl'] - m_Cl)).sum() < TOL
