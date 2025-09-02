def test_traj_bytesio():
    from ase_koopmans.io import iread, write
    from ase_koopmans.build import bulk
    from ase_koopmans.collections import g2
    import io

    images = [bulk('Si') + bulk('Fe')] + list(g2)

    buf = io.BytesIO()
    write(buf, images, format='traj')
    txt = buf.getvalue()

    buf = io.BytesIO()
    buf.write(txt)
    buf.seek(0)

    images2 = list(iread(buf, format='traj'))

    for atoms1, atoms2 in zip(images, images2):
        assert atoms1 == atoms2
