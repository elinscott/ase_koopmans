from ase.data import chemical_symbols, reference_states
from ase.build import bulk


for Z, ref in enumerate(reference_states):
    if ref is None:
        continue

    structure = ref['symmetry']
    if structure == 'diatom':
        continue  # not bulk
    if structure == 'atom':
        continue  # not bulk
    if structure == 'monoclinic':
        continue  # not implemented with bulk()
    if structure == 'cubic':
        continue  # (What is the meaning of the cubic (not sc) structures?)

    lat_map = dict(fcc='FCC',
                   bcc='BCC',
                   hcp='HEX',
                   tetragonal='TET',
                   diamond='FCC',
                   sc='CUB',
                   orthorhombic='ORC',
                   rhombohedral='RHL')

    sym = chemical_symbols[Z]
    atoms = bulk(sym)
    lat = atoms.cell.get_bravais_lattice()
    print(Z, atoms.symbols[0], structure, lat, atoms.cell.lengths())
    par1 = lat.tocell().niggli_reduce()[0].cellpar()
    par2 = atoms.cell.niggli_reduce()[0].cellpar()
    assert abs(par2 - par1).max() < 1e-10
    assert lat_map[structure] == lat.name

    if lat.name in ['RHL']:
        continue

    orth_atoms = bulk(sym, orthorhombic=True)
    orc_lat = orth_atoms.cell.get_bravais_lattice()
    angles = orc_lat.cellpar()[3:]
    assert abs(angles - 90).max() < 1e-10

    if lat.name in ['HEX', 'TET', 'ORC']:
        continue

    cub_atoms = bulk(sym, cubic=True)
    cub_lat = cub_atoms.cell.get_bravais_lattice()
    assert cub_lat.name == 'CUB', cub_lat
