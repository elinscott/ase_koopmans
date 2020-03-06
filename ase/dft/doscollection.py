from abc import ABCMeta
from collections.abc import Sequence as Sequence_abc
from functools import singledispatch
from typing import (Any, Iterable, List, Optional,
                    overload, Sequence, Tuple, Union)

import numpy as np
from ase.dft.dosdata import DOSData, RawDOSData, GridDOSData, Info


class DOSCollection(Sequence_abc, metaclass=ABCMeta):
    """Abstract base class for a collection of DOSData objects"""
    def __init__(self, dos_series: Iterable[DOSData]) -> None:
        self._data = list(dos_series)

    def sample(self,
               x: Sequence[float],
               width: float = 0.1,
               smearing: str = 'Gauss') -> np.ndarray:
        """Sample the DOS data at chosen points, with broadening

        This samples the underlying DOS data in the same way as the .sample()
        method of those DOSData items, returning a 2-D array with columns
        corresponding to x and rows corresponding to the collected data series.

        Args:
            x: energy values for sampling
            width: Width of broadening kernel
            smearing: selection of broadening kernel (only "Gauss" is currently
                supported)

        Returns:
            Weights sampled from a broadened DOS at values corresponding to x,
            in rows corresponding to DOSData entries contained in this object
        """
        return np.asarray([data.sample(x, width=width, smearing=smearing)
                           for data in self])

    def sample_grid(self,
                    npts: int,
                    xmin: float = None,
                    xmax: float = None,
                    padding: float = 3,
                    width: float = 0.1,
                    smearing: str = 'Gauss',
                    ) -> Tuple[Sequence[float], np.ndarray]:
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
            xmin = (min(min(data.get_energies()) for data in self)
                    - (padding * width))
        if xmax is None:
            xmax = (max(max(data.get_energies()) for data in self)
                    + (padding * width))
        x = np.linspace(xmin, xmax, npts)
        return x, self.sample(x, width=width, smearing=smearing)

    @classmethod
    def from_data(cls,
                  x: Sequence[float],
                  weights: Sequence[Sequence[float]],
                  info: Sequence[Info] = None) -> 'DOSCollection':
        """Create a DOSCollection from data sharing a common set of energies

        This is a convenience method to be used when all the DOS data in the
        collection has a common x-axis. There is no performance advantage in
        using this method for the generic DOSCollection, but for
        GridDOSCollection it is more efficient.

        Args:
            x: common set x-values for input data
            weights: array of DOS weights with rows corresponding to different
                datasets
            info: sequence of info dicts corresponding to weights rows.

        Returns:
            Collection of DOS data (in RawDOSData format)
        """

        info = cls._check_weights_and_info(weights, info)

        return cls(RawDOSData(x, row_weights, row_info)
                   for row_weights, row_info in zip(weights, info))

    @staticmethod
    def _check_weights_and_info(weights: Sequence[Sequence[float]],
                                info: Optional[Sequence[Info]],
                                ) -> Sequence[Info]:
        if info is None:
            info = [{}] * len(weights)
        else:
            if len(info) != len(weights):
                raise ValueError("Length of info must match number of rows in "
                                 "weights")
        return info

    @overload  # noqa F811
    def __getitem__(self, item: int) -> DOSData:
        ...

    @overload  # noqa F811
    def __getitem__(self, item: slice) -> List[DOSData]:
        ...

    def __getitem__(self, item): # noqa F811
        if isinstance(item, (int, slice)):
            return self._data[item]
        else:
            raise TypeError("index in DOSCollection must be an integer or "
                            "slice")

    def __len__(self) -> int:
        return len(self._data)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, type(self)):
            return False
        elif not len(self) == len(other):
            return False
        else:
            return all([a == b for a, b in zip(self, other)])

    def __add__(self,
                other: Union['DOSCollection', DOSData]) -> 'DOSCollection':
        """Join entries between two DOSCollection objects of the same type

        It is also possible to add a single DOSData object without wrapping it
        in a new collection: i.e.

          DOSCollection([dosdata1]) + DOSCollection([dosdata2])

        or

          DOSCollection([dosdata1]) + dosdata2

        will return

          DOSCollection([dosdata1, dosdata2])

        """
        return _add_to_collection(other, self)


@singledispatch
def _add_to_collection(other: DOSData,
                       collection: DOSCollection) -> DOSCollection:
    if isinstance(other, type(collection)):
        return type(collection)(list(collection) + list(other))
    else:
        raise TypeError("Only DOSCollection objects of the same type may "
                        "be joined with '+'.")


@_add_to_collection.register(DOSData)
def _(other: DOSData, collection: DOSCollection) -> DOSCollection:
    """Return a new DOSCollection with an additional DOSData item"""
    if isinstance(other, DOSData):
        return type(collection)(list(collection) + [other])


class RawDOSCollection(DOSCollection):
    def __init__(self, dos_series: Iterable[RawDOSData]) -> None:
        super().__init__(dos_series)
        for dos_data in self:
            if not isinstance(dos_data, RawDOSData):
                raise TypeError("RawDOSCollection can only store "
                                "RawDOSData objects.")


class GridDOSCollection(DOSCollection):
    def __init__(self, dos_series: Iterable[GridDOSData]) -> None:
        dos_list = list(dos_series)
        self._energies = dos_list[0].get_energies()
        self._weights = np.empty((len(dos_list), len(self._energies)), float)
        self._info = []

        for i, dos_data in enumerate(dos_list):
            if not isinstance(dos_data, GridDOSData):
                raise TypeError("GridDOSCollection can only store "
                                "GridDOSData objects.")
            if (dos_data.get_energies().shape != self._energies.shape
                or not np.allclose(dos_data.get_energies(), self._energies)):
                raise ValueError("All GridDOSData objects in GridDOSCollection"
                                 " must have the same energy axis.")
            self._weights[i, :] = dos_data.get_weights()
            self._info.append(dos_data.info)

    def __len__(self) -> int:
        return self._weights.shape[0]

    @overload  # noqa F811
    def __getitem__(self, item: int) -> DOSData:
        ...

    @overload  # noqa F811
    def __getitem__(self, item: slice) -> List[DOSData]:
        ...

    def __getitem__(self, item):  # noqa F811
        if isinstance(item, int):
            return GridDOSData(self._energies, self._weights[item, :],
                               info=self._info[item])
        elif isinstance(item, slice):
            return [self[i] for i in range(len(self))[item]]
        else:
            raise TypeError("index in DOSCollection must be an integer or "
                            "slice")

    @classmethod
    def from_data(cls,
                  x: Sequence[float],
                  weights: Sequence[Sequence[float]],
                  info: Sequence[Info] = None) -> 'DOSCollection':
        """Create a GridDOSCollection from data with a common set of energies

        This convenience method may also be more efficient as it limits
        redundant copying/checking of the data.

        Args:
            x: common set x-values for input data
            weights: array of DOS weights with rows corresponding to different
                datasets
            info: sequence of info dicts corresponding to weights rows.

        Returns:
            Collection of DOS data (in RawDOSData format)
        """

        weights_array = np.asarray(weights, dtype=float)
        if len(weights_array.shape) != 2:
            raise IndexError("Weights must be a 2-D array or nested sequence")
        if weights_array.shape[0] < 1:
            raise IndexError("Weights cannot be empty")
        if weights_array.shape[1] != len(x):
            raise IndexError("Length of weights rows must equal size of x")

        info = cls._check_weights_and_info(weights, info)

        dos_collection = cls([GridDOSData(x, weights_array[0])])
        dos_collection._weights = weights_array
        dos_collection._info = list(info)

        return dos_collection
