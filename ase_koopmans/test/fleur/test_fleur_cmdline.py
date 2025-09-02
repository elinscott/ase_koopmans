def test_fleur_cmdline(cli):
    cli.shell('ase_koopmans build -x fcc -a 4.04 Al | ase_koopmans run fleur -p kpts=3.0,xc=PBE',
              'fleur')
