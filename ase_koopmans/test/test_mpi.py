def test_mpi():
    """Try to import all ASE modules and check that ase_koopmans.parallel.world has not
    been used.  We want to delay use of world until after MPI4PY has been
    imported.

    We run the test in a subprocess so that we have a clean Python interpreter."""

    import sys
    from subprocess import run

    # Should cover most of ASE:
    modules = ['ase_koopmans.optimize',
               'ase_koopmans.db',
               'ase_koopmans.gui']

    imports = 'import ' + ', '.join(modules)

    run([sys.executable,
         '-c',
         '{imports}; from ase_koopmans.parallel import world; assert world.comm is None'
         .format(imports=imports)],
        check=True)
