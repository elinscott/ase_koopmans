def test_exciting():
    from ase_koopmans import Atoms
    from ase_koopmans.io import read, write
    from ase_koopmans.calculators.exciting import Exciting


    a = Atoms('N3O',
              [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0.5, 0.5, 0.5)],
              pbc=True)

    write('input.xml', a)
    b = read('input.xml')

    print(a)
    print(a.get_positions())
    print(b)
    print(b.get_positions())

    Exciting(dir='excitingtestfiles',
             kpts=(4, 4, 3),
             # bin='/fshome/chm/git/exciting/bin/excitingser',
             maxscl=3)
    # maybe do something???
