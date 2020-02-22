# Refactor of DOS-like data objects
# towards replacing ase.dft.dos and ase.dft.pdos
from abc import ABCMeta, abstractmethod
import logging
from typing import Dict, Sequence, Tuple

from matplotlib.axes import Axes
import numpy as np
from scipy.signal import convolve


# For now we will be strict about Info and say it has to be str->str. Perhaps
# later we will allow other types that have reliable comparison operations.
Info = Dict[str, str]


class DOSData(metaclass=ABCMeta):
    """Abstract base class for a single series of DOS-like data"""
    def __init__(self,
                 info: Info = None) -> None:
        if info is None:
            self.info = {}
        elif isinstance(info, dict):
            self.info = info
        else:
            raise TypeError("Info must be a dict or None")

    @abstractmethod
    def get_energies(self) -> Sequence[float]:
        """Get energy data stored in this object"""
        raise NotImplementedError

    @abstractmethod
    def get_weights(self) -> Sequence[float]:
        """Get DOS weights stored in this object"""
        raise NotImplementedError

    def sample(self,
               x: Sequence[float],
               width: float = 0.1,
               smearing: str = 'Gauss') -> np.ndarray:
        """Sample the DOS data at chosen points, with broadening"""

        self._check_positive_width(width)
        weights_grid = np.dot(
            self.get_weights(),
            self._delta(np.asarray(x),
                        np.asarray(self.get_energies())[:, np.newaxis],
                        width,
                        smearing=smearing))
        return weights_grid

    @staticmethod
    def _delta(x: np.ndarray,
               x0: np.ndarray,
               width: float,
               smearing: str = 'Gauss') -> Sequence[Sequence[float]]:
        """Return a delta-function centered at 'x0'.

        This function is used with numpy broadcasting; if x is a row and x0 is
        a column vector, the returned data will be a 2D array with each row
        corresponding to a different delta center.
        """
        if smearing.lower() == 'gauss':
            x1 = -0.5 * ((x - x0) / width)**2
            return np.exp(x1) / (np.sqrt(2 * np.pi) * width)
        else:
            msg = 'Requested smearing type not recognized. Got {}'.format(
                smearing)
            raise ValueError(msg)

    @staticmethod
    def _check_positive_width(width):
        if width <= 0.0:
            msg = 'Cannot add 0 or negative width smearing'
            raise ValueError(msg)

    def sample_grid(self,
                    npts: int,
                    xmin: float = None,
                    xmax: float = None,
                    padding: float = 3,
                    width: float = 0.1,
                    smearing: str = 'Gauss',
                    ) -> Tuple[Sequence[float], Sequence[float]]:
        """Sample the DOS data on an evenly-spaced energy grid

        Args:
            npts: Number of sampled points
            xmin: Minimum sampled x value; if unspecified, a default is chosen
            xmax: Maximum sampled x value; if unspecified, a default is chosen
            padding: If xmin/xmax is unspecified, default value will be padded
                by padding * width to avoid cutting off peaks.
            width: Width of broadening kernel, passed to self.sample()
            smearing: selection of broadening kernel, passed to self.sample()

        Returns:
            (x-values, sampled DOS)
        """

        if xmin is None:
            xmin = min(self.get_energies()) - (padding * width)
        if xmax is None:
            xmax = max(self.get_energies()) + (padding * width)
        x = np.linspace(xmin, xmax, npts)
        return x, self.sample(x, width=width, smearing=smearing)

    def plot_dos(self,
                 npts: int = 1000,
                 xmin: float = None,
                 xmax: float = None,
                 width: float = 0.1,
                 smearing: str = 'Gauss',
                 ax: Axes = None,
                 show: bool = False,
                 filename: str = None,
                 mplargs: dict = None) -> Axes:
        """Simple 1-D plot of DOS data, resampled onto a grid

        If the special key 'label' is present in self.info, this will be set
        as the label for the plotted line (unless overruled in mplargs). The
        label is only seen if a legend is added to the plot (i.e. by calling
        `ax.legend()`).

        Args:
            npts, xmin, xmax: output data range, as passed to self.sample_grid
            width: Width of broadening kernel, passed to self.sample()
            smearing: selection of broadening kernel, passed to self.sample()
            ax: existing Matplotlib axes object. If not provided, a new figure
                with one set of axes will be created using Pyplot
            show: show the figure on-screen
            filename: if a path is given, save the figure to this file
            mplargs: additional arguments to pass to matplotlib plot command
                (e.g. {'linewidth': 2} for a thicker line).


        Returns:
            Plotting axes. If "ax" was set, this is the same object.
        """

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        if mplargs is None:
            mplargs = {}

        x, y = self.sample_grid(npts, xmin=xmin, xmax=xmax,
                                width=width, smearing=smearing)
        ax.plot(x, y, **mplargs)

        if show:
            fig.show()
        if filename is not None:
            fig.savefig(filename)

        return ax


class RawDOSData(DOSData):
    """A collection of weighted delta functions which sum to form a DOS

    This is an appropriate data container for density-of-states (DOS) or
    spectral data where the energy data values not form a known regular
    grid. The data may be plotted or resampled for further analysis using the
    sample(), sample_grid() and plot() methods. Multiple weights at the same
    energy value will _only_ be combined in output data, and data stored in
    RawDOSData is never resampled. A plot_deltas() function is also provided
    which plots the raw data.

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

    def get_energies(self) -> np.ndarray:
        return self._data[0, :].copy()

    def get_weights(self) -> np.ndarray:
        return self._data[1, :].copy()

    def __add__(self, other: 'RawDOSData') -> 'RawDOSData':
        if not isinstance(other, RawDOSData):
            raise TypeError("RawDOSData can only be combined with other "
                            "RawDOSData objects")

        # Take intersection of metadata (i.e. only common entries are retained)
        new_info = dict(set(self.info.items()) & set(other.info.items()))

        # Concatenate the energy/weight data
        new_data = np.concatenate((self._data, other._data), axis=1)

        new_object = RawDOSData([], [], info=new_info)
        new_object._data = new_data

        return new_object

    def plot_deltas(self,
                    ax: Axes = None,
                    show: bool = False,
                    filename: str = None,
                    mplargs: dict = None) -> Axes:
        """Simple plot of sparse DOS data as a set of delta functions

        Items at the same x-value can overlap and will not be summed together

        Args:
            ax: existing Matplotlib axes object. If not provided, a new figure
                with one set of axes will be created using Pyplot
            show: show the figure on-screen
            filename: if a path is given, save the figure to this file
            mplargs: additional arguments to pass to matplotlib Axes.vlines
                command (e.g. {'linewidth': 2} for a thicker line).

        Returns:
            Plotting axes. If "ax" was set, this is the same object.
        """
        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        if mplargs is None:
            mplargs = {}

        ax.vlines(self.get_energies(), 0, self.get_weights(), **mplargs)

        if show:
            fig.show()
        if filename is not None:
            fig.savefig(filename)

        return ax


class GridDOSData(DOSData):
    """A collection of regularly-sampled data which represents a DOS

    This is an appropriate data container for density-of-states (DOS) or
    spectral data where the intensity values form a regular grid. This is
    generally the result of sampling or integrating into discrete bins,
    rather than a collection of unique states. The data may be plotted or
    resampled for further analysis using the sample(), sample_grid() and plot()
    methods.

    Metadata may be stored in the info dict, in which keys and values must be
    strings. This data is used for selecting and combining multiple DOSData
    objects in a DOSCollection object.

    When RawDOSData objects are combined with the addition operator::

      big_dos = raw_dos_1 + raw_dos_2

    the weights data is _summed_ (requiring a consistent energy grid) and the
    new info dictionary consists of the _intersection_ of the inputs: only
    key-value pairs that were common to both of the input objects will be
    retained in the new combined object. For example:

    (GridDOSData([0.1, 0.2, 0.3], [y1, y2, y3],
                 info={'symbol': 'O', 'index': '1'})
     + GridDOSData([0.1, 0.2, 0.3], [y4, y5, y6],
                   info={'symbol': 'O', 'index': '2'})

    will yield the equivalent of

    GridDOSData([0.1, 0.2, 0.3], [y1+y4, y2+y5, y3+y6], info={'symbol': 'O'})

    """
    def __init__(self,
                 energies: Sequence[float],
                 weights: Sequence[float],
                 info: Info = None) -> None:
        super().__init__(info=info)

        n_entries = len(energies)
        if not np.allclose(energies,
                           np.linspace(energies[0], energies[-1], n_entries)):
            raise ValueError("Energies must be an evenly-spaced 1-D grid")

        if len(weights) != n_entries:
            raise ValueError("Energies and weights must be the same length")

        # Internally store the data as a np array with two rows; energy, weight
        self._data = np.empty((2, n_entries), dtype=float, order='C')
        self._data[0, :] = energies
        self._data[1, :] = weights

        self.sigma_cutoff = 3

    def get_energies(self) -> np.ndarray:
        return self._data[0, :].copy()

    def get_weights(self) -> np.ndarray:
        return self._data[1, :].copy()

    def _check_spacing(self, width):
        current_spacing = self._data[0, 1] - self._data[0, 0]
        if width < (2 * current_spacing):
            logging.warning(
                "The broadening width is small compared to the original "
                "sampling density. The results are unlikely to be smooth.")

    def sample(self,
               x: Sequence[float],
               width: float = 0.1,
               smearing: str = 'Gauss') -> np.ndarray:
        self._check_spacing(width)
        return super().sample(x=x, width=width, smearing=smearing)

    def sample_grid(self,
                    npts: int,
                    xmin: float = None,
                    xmax: float = None,
                    padding: float = 3,
                    width: float = 0.1,
                    smearing: str = 'Gauss'
                    ) -> Tuple[Sequence[float], Sequence[float]]:
        """Resample DOS data, broadening with convolution

        Args:
            npts: Number of sampled points
            xmin: Minimum sampled x value; if unspecified, a default is chosen
            xmax: Maximum sampled x value; if unspecified, a default is chosen
            padding: If xmin/xmax is unspecified, default value will be padded
                by padding * width to avoid cutting off peaks.
            width: Width of broadening kernel, passed to self.sample()
            smearing: selection of broadening kernel, passed to self.sample()

        Returns:
            (x-values, sampled DOS)

        The broadening kernel is closely approximated within a finite region:
        The function is truncated at a distance of 3*width; for
        a Gaussian this corresponds to a value of ~10^-4. The kernel is shifted
        and rescaled to pass through zero at this point, maintaining unity
        intensity.

        """
        self._check_positive_width(width)
        self._check_spacing(width)

        energies = self.get_energies()
        if xmin is None:
            xmin = min(energies) - (padding * width)
        if xmax is None:
            xmax = max(energies) + (padding * width)

        new_energies = np.linspace(xmin, xmax, npts)
        spacing = new_energies[1] - new_energies[0]
        # Define histogram bins as evenly surrounding x points
        bins = np.linspace(xmin - spacing / 2, xmax + spacing / 2, npts + 1)

        kernel_range_in_pts = int(np.ceil(self.sigma_cutoff * width / spacing))
        kernel_range = kernel_range_in_pts * spacing
        kernel_x_values = np.linspace(-kernel_range, kernel_range,
                                      kernel_range_in_pts * 2 + 1)

        if smearing.lower() == 'gauss':
            # Normalisation term of Gaussian is removed as we will next
            # normalise "the hard way" to account for truncation.
            kernel = np.exp(-0.5 * (kernel_x_values / width)**2)

        else:
            msg = 'Requested smearing type not recognized. Got {}'.format(
                smearing)
            raise ValueError(msg)

        # Remove "step" at cutoff
        kernel -= kernel[0]
        # Normalise such that peaks sum to 1
        kernel /= kernel.sum()
        
        # Re-bin the data to the new x-grid, rescaling for consistent density
        input_spacing = energies[1] - energies[0]
        weight_scale = input_spacing / spacing
        weights, _ = np.histogram(energies, bins=bins,
                                  weights=self.get_weights() * weight_scale,
                                  density=False)

        # Convolve with kernel on same grid
        return new_energies, convolve(weights, kernel, mode='same')

    def __add__(self, other: 'GridDOSData') -> 'GridDOSData':
        # This method uses direct access to the mutable energy and weights data
        # (self._data) to avoid redundant copying operations. The __init__
        # method of GridDOSData will write this to a new array, so on this
        # occasion it is safe to pass references to the mutable data.

        if not isinstance(other, GridDOSData):
            raise TypeError("GridDOSData can only be combined with other "
                            "GridDOSData objects")
        if not np.allclose(self._data[0, :], other.get_energies()):
            raise ValueError("Cannot add GridDOSData objects with different "
                             "energy grids.")

        # Take intersection of metadata (i.e. only common entries are retained)
        new_info = dict(set(self.info.items()) & set(other.info.items()))

        # Concatenate the energy/weight data
        new_weights = self._data[1, :] + other.get_weights()

        new_object = GridDOSData(self._data[0, :], new_weights,
                                 info=new_info)
        return new_object
