def test_ag(cli):
    from ase_koopmans import Atoms
    from ase_koopmans.io import write

    write('x.json', Atoms('X'))

    # Make sure ASE's gui can run in terminal mode without $DISPLAY and tkinter:
    cli.shell('ase_koopmans -T gui --terminal -n "id=1" x.json')
