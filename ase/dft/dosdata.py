# Refactor of DOS-like data objects
# towards replacing ase.dft.dos and ase.dft.pdos

from matplotlib.axes import Axes
import numpy as np
from typing import Dict, Sequence
Info = Dict[str, str]  # For now we will be strict about this. Perhaps later we
                       # will allow other types that have reliable comparison
                       # operations.

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
