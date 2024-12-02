def test_aims_cmdline(cli):
    # warning! parameters are not converged - only an illustration!
    cli.shell("""ase_koopmans build -x bcc -a 3.6 Li | \
    ase_koopmans run aims -s 0.3 -p \
    kpts=1.5,xc=LDA,sc_accuracy_rho=5.e-2,relativistic=none,compute_analytical_stress=True,sc_accuracy_forces=5.e-1""", 'aims')
