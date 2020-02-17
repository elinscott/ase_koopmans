# Refactor of DOS-like data objects
# towards replacing ase.dft.dos and ase.dft.pdos

from matplotlib.axes import Axes
import numpy as np
from typing import Dict, Sequence
Info = Dict[str, str]  # For now we will be strict about this. Perhaps later we
                       # will allow other types that have reliable comparison
                       # operations.

def RawDOSData(DOSData):
    """A collection of weighted delta functions which sum to form a DOS

    This is an appropriate data container for DOS or spectral data where the
    energy data values not form a known regular grid. The data may be plotted
    or resampled for further analysis using the sample(), sample_grid() and
    plot() methods. Multiple weights at the same energy value will _only_ be
    combined in output data, and data stored in RawDOSData is never resampled.

    Metadata may be stored in the info dict, in which keys and values must be
    strings. This data is used for selecting and combining multiple DOSData
    objects in a DOSCollection object.

    When RawDOSData objects are combined with the addition operator::

      big_dos = raw_dos_1 + raw_dos_2

    the energy and weights data is _concatenated_ (i.e. combined without
    sorting or replacement) and the new info dictionary consists of the
    _intersection_ of the inputs: only key-value pairs that were common to both
    of the input objects will be retained in the new combined object. For 
    example:

    (RawDOSData([x1], [y1], info={'symbol': 'O', 'index': '1'})
     + RawDOSData([x2], [y2], info={'symbol': 'O', 'index': '2'})

    will yield the equivalent of

    RawDOSData([x1, x2], [y1, y2], info={'symbol': 'O'})

    """
    def __init__(self,
                 energies: Sequence[float],
                 weights: Sequence[float],
                 info: Info = None) -> None:
        super().__init__(info=info)

        n_entries = len(energies)
        if len(weights) != n_entries:
            raise ValueError("Energies and weights must be the same length")

        # Internally store the data as a np array with two rows; energy, weight
        self._data = np.empty((2, n_entries), dtype=float, order='C')
        self._data[0, :] = energies
        self._data[1, :] = weights

    def get_energies(self) -> Sequence[float]:
        return self._data[0, :].copy()

    def get_weights(self) -> Sequence[float]:
        return self._data[1, :].copy()

    def sample(self, x, width=0.1, smearing='Gauss'):
        if width <= 0.0:
            msg = 'Cannot add 0 or negative width smearing'
            raise ValueError(msg)

        weights_grid = np.dot(self.get_weights(),
                              self.delta(x,
                                         self.get_energies()[:, np.newaxis],
                                         width,
                                         smearing=smearing))
        return weights_grid    

    @staticmethod
    def _delta(x, x0, width, smearing='Gauss'):
        """Return a delta-function centered at 'x0'."""
        if smearing.lower() == 'gauss':
            x1 = -((x - x0) / width)**2
            return np.exp(x1) / (np.sqrt(np.pi) * width)
        else:
            msg = 'Requested smearing type not recognized. Got {}'.format(
                smearing)
            raise ValueError(msg)



def DOSData(object):
    """Abstract base class for a single series of DOS-like data"""
    def __init__(self,
                 info: Info = None) -> None:
        if info is None:
            self.info = {}
        else:
            self.info = info

    def get_energies(self) -> sequence[float]:
        """Get energy data stored in this object"""
        raise NotImplementedError

    def get_weights(self) -> sequence[float]:
        """Get DOS weights stored in this object"""
        raise NotImplementedError

    def sample(self, x: Sequence(float), **broadening_args) -> sequence[float]:
        """Sample the DOS data at chosen points, with broadening"""
        raise NotImplementedError

    def sample_grid(self,
                    npts: int
                    xmin: float = None,
                    xmax: float = None,
                    **broadening_args
                    ) -> Tuple[sequence[float], sequence[float]]:
        """Sample the DOS data on an evenly-spaced energy grid

        Args:
            npts: Number of sampled points
            xmin: Minimum sampled x value; if unspecified, a default is chosen
            xmax: Maximum sampled x value; if unspecified, a default is chosen
            **broadening_args: broadening options passed to self.sample()

        Returns:
            (x-values, sampled DOS)
        """

        if xmin is None:
            xmin = min(self.get_energies())
        if xmax is None:
            xmax = max(self.get_weights())
        x = np.linspace(xmin, xmax, npts)
        return x, self.sample(x, **broadening_args)

    def plot_dos(self,
                 npts: int = 1000,
                 xmin: float = None,
                 xmax: float = None
                 ax: Axes = None,
                 show: bool = False,
                 filename: str = None,
                 mplargs: dict = None,
                 **broadening_args) -> Axes:
        """Simple 1-D plot of DOS data, resampled onto a grid

        If the special key 'label' is present in self.info, this will be set
        as the label for the plotted line (unless overruled in mplargs). The
        label is only seen if a legend is added to the plot (i.e. by calling
        `ax.legend()`).

        Args:
            npts, xmin, xmax: output data range, as passed to self.sample_grid
            ax: existing Matplotlib axes object. If not provided, a new figure
                with one set of axes will be created using Pyplot
            show: show the figure on-screen
            filename: if a path is given, save the figure to this file
            mplargs: additional arguments to pass to matplotlib plot command 
                (e.g. {'linewidth': 2} for a thicker line.
            **broadening_args: remaining keyword arguments are used for
                  sampling/broadening and should be supported by self.sample()

        Returns:
            Plotting axes. If "ax" was set, this is the same object.
        """

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plot.subplots()
        else:
            fig = ax.get_figure()

        x, y = self.sample_grid(npts, xmin=xmin, xmax=xmax, **broadening_args)
        ax.plot(x, y, **mplargs)

        if show:
            fig.show()
        if filename is not None:
            fig.savefig(filename)

        return ax
