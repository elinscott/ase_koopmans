def test_ag():
    from ase import Atoms
    from ase.io import write
    from ase.test import cli


    write('x.json', Atoms('X'))

    # Make sure ASE's gui can run in terminal mode without $DISPLAY and tkinter:
    cli('ase -T gui --terminal -n "id=1" x.json')
