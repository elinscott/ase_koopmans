.. module:: ase.spectrum.dosdata
   :synopsis: Density of states data objects

======================
Density of states data
======================
*The dosdata and doscollection modules contained in ase.spectrum form*
*a new framework, intended to replace the DOS implementations in*
*ase.dft.dos and ase.dft.pdos.*

This module provides the RawDOSData and GridDOSData objects that
contain a single series of spectral contributions with corresponding
energies. This can then be resampled with broadening on-the-fly to
create typical density-of-states plots.

Any DOS data object includes a ``.info`` dictionary of metadata. The
classes in :obj:`ase.spectrum.doscollection` can store multiple
DOSData objects and have methods for selecting and merging data series
based on this metadata.

DOSData is an abstract base class and should not be used
directly. RawDOSData is intended for "sparse" data such as an exact
set of eigenvalues for electronic states. Precision errors related to
binning/broadening are avoided by deferring binning until data is
plotted or extracted using ``get_dos()``. GridDOSData stores data on a
regular energy series, and is suitable for data that was imported in
this format (e.g. a Bloechl-corrected tetrahedron-integrated DOS on a
dense k-point mesh).

More details
------------

.. automodule:: ase.spectrum.dosdata
   :members:
   :noindex:
