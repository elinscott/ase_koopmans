from ase.lattice import RHL


def test_ase_bandstructure(cli):
    lat = RHL(3., 70.0)
    path = lat.bandpath()
    bs = path.free_electron_band_structure()
    bs.write('bs.json')

    cli.ase('band-structure bs.json --output bs.png')
    # If the CLI tool gave a text output, we could verify it.


# Note: We don't have proper testing of --points, --range etc.  We
# test here on JSON input but the tool is in principle supposed to
# work on other formats, too (only gpw though as of now though).
