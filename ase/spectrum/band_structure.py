import numpy as np

import ase  # Annotations
from ase.utils import jsonable
from ase.calculators.calculator import PropertyNotImplementedError


def calculate_band_structure(atoms, path=None, scf_kwargs=None,
                             bs_kwargs=None, kpts_tol=1e-6, cell_tol=1e-6):
    """Calculate band structure.

    The purpose of this function is to abstract a band structure calculation
    so the workflow does not depend on the calculator.

    First trigger SCF calculation if necessary, then set arguments
    on the calculator for band structure calculation, then return
    calculated band structure.

    The difference from get_band_structure() is that the latter
    expects the calculation to already have been done."""
    if path is None:
        path = atoms.cell.bandpath()

    from ase.lattice import celldiff  # Should this be a method on cell?
    if any(path.cell.any(1) != atoms.pbc):
        raise ValueError('The band path\'s cell, {}, does not match the '
                         'periodicity {} of the atoms'
                         .format(path.cell, atoms.pbc))
    cell_err = celldiff(path.cell, atoms.cell.uncomplete(atoms.pbc))
    if cell_err > cell_tol:
        raise ValueError('Atoms and band path have different unit cells.  '
                         'Please reduce atoms to standard form.  '
                         'Cell lengths and angles are {} vs {}'
                         .format(atoms.cell.cellpar(), path.cell.cellpar()))

    calc = atoms.calc
    if calc is None:
        raise ValueError('Atoms have no calculator')

    if scf_kwargs is not None:
        calc.set(**scf_kwargs)

    # Proposed standard mechanism for calculators to advertise that they
    # use the bandpath keyword to handle band structures rather than
    # a double (SCF + BS) run.
    use_bandpath_kw = getattr(calc, 'accepts_bandpath_keyword', False)
    if use_bandpath_kw:
        calc.set(bandpath=path)
        atoms.get_potential_energy()
        return calc.band_structure()

    atoms.get_potential_energy()

    if hasattr(calc, 'get_fermi_level'):
        # What is the protocol for a calculator to tell whether
        # it has fermi_energy?
        eref = calc.get_fermi_level()
    else:
        eref = 0.0

    if bs_kwargs is None:
        bs_kwargs = {}

    calc.set(kpts=path, **bs_kwargs)
    calc.results.clear()  # XXX get rid of me

    # Calculators are too inconsistent here:
    # * atoms.get_potential_energy() will fail when total energy is
    #   not in results after BS calculation (Espresso)
    # * calc.calculate(atoms) doesn't ask for any quantity, so some
    #   calculators may not calculate anything at all
    # * 'bandstructure' is not a recognized property we can ask for
    try:
        atoms.get_potential_energy()
    except PropertyNotImplementedError:
        pass

    ibzkpts = calc.get_ibz_k_points()
    kpts_err = np.abs(path.kpts - ibzkpts).max()
    if kpts_err > kpts_tol:
        raise RuntimeError('Kpoints of calculator differ from those '
                           'of the band path we just used; '
                           'err={} > tol={}'.format(kpts_err, kpts_tol))

    bs = get_band_structure(atoms, path=path, reference=eref)
    return bs


def get_band_structure(atoms=None, calc=None, path=None, reference=None):
    """Create band structure object from Atoms or calculator."""
    # path and reference are used internally at the moment, but
    # the exact implementation will probably change.  WIP.
    #
    # XXX We throw away info about the bandpath when we create the calculator.
    # If we have kept the bandpath, we can provide it as an argument here.
    # It would be wise to check that the bandpath kpoints are the same as
    # those stored in the calculator.
    atoms = atoms if atoms is not None else calc.atoms
    calc = calc if calc is not None else atoms.calc

    kpts = calc.get_ibz_k_points()

    energies = []
    for s in range(calc.get_number_of_spins()):
        energies.append([calc.get_eigenvalues(kpt=k, spin=s)
                         for k in range(len(kpts))])
    energies = np.array(energies)

    if path is None:
        from ase.dft.kpoints import resolve_custom_points, find_bandpath_kinks
        path = atoms.cell.bandpath(npoints=0)
        # Kpoints are already evaluated, we just need to put them into
        # the path (whether they fit our idea of what the path is, or not).
        #
        # Depending on how the path was established, the kpoints might
        # be valid high-symmetry points, but since there are multiple
        # high-symmetry points of each type, they may not coincide
        # with ours if the bandpath was generated by another code.
        #
        # Here we hack it so the BandPath has proper points even if they
        # come from some weird source.
        #
        # This operation (manually hacking the bandpath) is liable to break.
        # TODO: Make it available as a proper (documented) bandpath method.
        kinks = find_bandpath_kinks(atoms.cell, kpts, eps=1e-5)
        pathspec = resolve_custom_points(kpts[kinks], path.special_points, eps=1e-5)
        path._kpts = kpts
        path._path = pathspec

    # XXX If we *did* get the path, now would be a good time to check
    # that it matches the cell!  Although the path can only be passed
    # because we internally want to not re-evaluate the Bravais
    # lattice type.  (We actually need an eps parameter, too.)

    if reference is None:
        # Fermi level should come from the GS calculation, not the BS one!
        reference = calc.get_fermi_level()

    if reference is None:
        # Fermi level may not be available, e.g., with non-Fermi smearing.
        # XXX Actually get_fermi_level() should raise an error when Fermi
        # level wasn't available, so we should fix that.
        reference = 0.0

    return BandStructure(path=path,
                         energies=energies,
                         reference=reference)


class BandStructurePlot:
    def __init__(self, bs):
        self.bs = bs
        self.ax = None
        self.xcoords = None
        self.show_legend = False

    def plot(self, ax=None, spin=None, emin=-10, emax=5, filename=None,
             show=False, ylabel=None, colors=None, label=None,
             spin_labels=['spin up', 'spin down'], loc=None, **plotkwargs):
        """Plot band-structure.

        spin: int or None
            Spin channel.  Default behaviour is to plot both spin up and down
            for spin-polarized calculations.
        emin,emax: float
            Maximum energy above reference.
        filename: str
            Write image to a file.
        ax: Axes
            MatPlotLib Axes object.  Will be created if not supplied.
        show: bool
            Show the image.
        """

        if self.ax is None:
            ax = self.prepare_plot(ax, emin, emax, ylabel)

        if spin is None:
            e_skn = self.bs.energies
        else:
            e_skn = self.bs.energies[spin, np.newaxis]

        if colors is None:
            if len(e_skn) == 1:
                colors = 'g'
            else:
                colors = 'yb'

        nspins = len(e_skn)

        for spin, e_kn in enumerate(e_skn):
            color = colors[spin]
            kwargs = dict(color=color)
            kwargs.update(plotkwargs)
            if nspins == 2:
                if label:
                    lbl = label + ' ' + spin_labels[spin]
                else:
                    lbl = spin_labels[spin]
            else:
                lbl = label
            ax.plot(self.xcoords, e_kn[:, 0], label=lbl, **kwargs)

            for e_k in e_kn.T[1:]:
                ax.plot(self.xcoords, e_k, **kwargs)

        self.show_legend = label is not None or nspins == 2
        self.finish_plot(filename, show, loc)

        return ax

    def plot_with_colors(self, ax=None, emin=-10, emax=5, filename=None,
                         show=False, energies=None, colors=None,
                         ylabel=None, clabel='$s_z$', cmin=-1.0, cmax=1.0,
                         sortcolors=False, loc=None, s=2):
        """Plot band-structure with colors."""

        import matplotlib.pyplot as plt

        if self.ax is None:
            ax = self.prepare_plot(ax, emin, emax, ylabel)

        shape = energies.shape
        xcoords = np.vstack([self.xcoords] * shape[1])
        if sortcolors:
            perm = colors.argsort(axis=None)
            energies = energies.ravel()[perm].reshape(shape)
            colors = colors.ravel()[perm].reshape(shape)
            xcoords = xcoords.ravel()[perm].reshape(shape)

        for e_k, c_k, x_k in zip(energies, colors, xcoords):
            things = ax.scatter(x_k, e_k, c=c_k, s=s,
                                vmin=cmin, vmax=cmax)

        cbar = plt.colorbar(things)
        cbar.set_label(clabel)

        self.finish_plot(filename, show, loc)

        return ax

    def prepare_plot(self, ax=None, emin=-10, emax=5, ylabel=None):
        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.figure().add_subplot(111)

        def pretty(kpt):
            if kpt == 'G':
                kpt = r'$\Gamma$'
            elif len(kpt) == 2:
                kpt = kpt[0] + '$_' + kpt[1] + '$'
            return kpt

        self.xcoords, label_xcoords, orig_labels = self.bs.get_labels()
        label_xcoords = list(label_xcoords)
        labels = [pretty(name) for name in orig_labels]

        i = 1
        while i < len(labels):
            if label_xcoords[i - 1] == label_xcoords[i]:
                labels[i - 1] = labels[i - 1] + ',' + labels[i]
                labels.pop(i)
                label_xcoords.pop(i)
            else:
                i += 1

        for x in label_xcoords[1:-1]:
            ax.axvline(x, color='0.5')

        ylabel = ylabel if ylabel is not None else 'energies [eV]'

        ax.set_xticks(label_xcoords)
        ax.set_xticklabels(labels)
        ax.set_ylabel(ylabel)
        ax.axhline(self.bs.reference, color='k', ls=':')
        ax.axis(xmin=0, xmax=self.xcoords[-1], ymin=emin, ymax=emax)
        self.ax = ax
        return ax

    def finish_plot(self, filename, show, loc):
        import matplotlib.pyplot as plt

        if self.show_legend:
            leg = self.ax.legend(loc=loc)
            leg.get_frame().set_alpha(1)

        if filename:
            plt.savefig(filename)

        if show:
            plt.show()


@jsonable('bandstructure')
class BandStructure:
    """A band structure consists of an array of eigenvalues and a bandpath.

    BandStructure objects support JSON I/O.
    """

    def __init__(self, path, energies, reference=0.0):
        self._path = path
        self._energies = np.asarray(energies)
        assert self.energies.shape[0] in [1, 2]  # spins x kpts x bands
        assert self.energies.shape[1] == len(path.kpts)
        assert np.isscalar(reference)
        self._reference = reference

    @property
    def energies(self) -> np.ndarray:
        """The energies of this band structure.

        This is a numpy array of shape (nspins, nkpoints, nbands)."""
        return self._energies

    @property
    def path(self) -> 'ase.dft.kpoints.BandPath':
        """The :class:`~ase.dft.kpoints.BandPath` of this band structure."""
        return self._path

    @property
    def reference(self) -> float:
        """The reference energy.

        Semantics may vary; typically a Fermi energy or zero,
        depending on how the band structure was created."""
        return self._reference

    def subtract_reference(self, reference=None) -> 'BandStructure':
        """Return new band structure with reference energy subtracted."""
        reference = reference if reference is not None else self.reference
        return BandStructure(self.path, self.energies - reference,
                             reference=0.0)

    def todict(self):
        return dict(path=self.path,
                    energies=self.energies,
                    reference=self.reference)

    def get_labels(self, eps=1e-5):
        """"See :func:`ase.dft.kpoints.labels_from_kpts`."""
        return self.path.get_linear_kpoint_axis(eps=eps)

    def plot(self, *args, **kwargs):
        """Plot this band structure."""
        bsp = BandStructurePlot(self)
        return bsp.plot(*args, **kwargs)

    def __repr__(self):
        return ('{}(path={!r}, energies=[{} values], reference={})'
                .format(self.__class__.__name__, self.path,
                        '{}x{}x{}'.format(*self.energies.shape),
                        self.reference))

    def __eq__(self, other):
        if isinstance(other, BandStructure):
            if not abs(self.reference - other.reference) < 1e-10:
                return False
            elif self.path != other.path:
                return False
            else:
                return np.allclose(self.energies, other.energies, atol=1e-10)
        return False
