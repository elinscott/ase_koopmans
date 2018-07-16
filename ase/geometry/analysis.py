"""Tools for analyzing instances of :class:`~ase.Atoms`
"""

import numpy as np
from ase.neighborlist import buildNeighborList, get_distance_matrix, get_distance_indices
from ase.ga.utilities import get_rdf
from ase import Atoms

#memory-friendly iterator based zip for python2
try:
    from itertools import izip as zip
except ImportError:
    pass

class Analysis(object):
    """Initialize Analysis class

    - ``images``: :class:`~ase.Atoms` object or list of such
    - ``nl``: None, :class:`~ase.neighborlist.NeighborList` object or list of such
    - ``**kwargs``: Arguments for constructing :class:`~ase.neighborlist.NeighborList` object if ``nl`` is None.

    The choice of ``bothways=True`` for the :class:`~ase.neighborlist.NeighborList` object
    will not influence the amount of bonds/angles/dihedrals you get, all are reported
    in both directions. Use the *unique*-labeled properties to get lists without
    duplicates.
    """

    def __init__(self, images, nl=None, **kwargs):
        self.images = images

        if isinstance(nl, list):
            assert len(nl) == self.nImages
            self._nl = nl
        elif nl is not None:
            self._nl = [ nl ]
        else:
            self._nl = [ buildNeighborList(self.images[0], **kwargs) ]

        self._adjacencyMatrix = None
        self._distanceMatrix = None
        self._allBonds = None
        self._allAngles = None
        self._allDihedrals = None

    def _delete_derived_info(self):
        del self.adjacency_matrix
        del self.distance_matrix
        del self.all_bonds
        del self.all_angles
        del self.all_dihedrals

    def _get_slice(self, imageIdx):
        """Return a slice from user input.
        Using *imageIdx* (can be integer or slice) the analyzed frames can be specified.
        If *imageIdx* is None, all frames will be analyzed.
        """
        #get slice from imageIdx
        if isinstance(imageIdx, int):
            sl = slice(imageIdx, imageIdx+1)
        elif isinstance(imageIdx, slice):
            sl = imageIdx
        elif imageIdx is None:
            sl = slice(0, None)
        else:
            raise ValueError("Unsupported type for imageIdx in ase.geometry.analysis.Analysis._get_slice")
        return sl

    @property
    def images(self):
        """Get images"""
        return self._images

    @images.setter
    def images(self, images):
        """Set images"""
        if isinstance(images, list):
            self._images = images
        else:
            self._images = [ images ]

        self._nImages = len(images)

    @images.deleter
    def images(self):
        """Delete images"""
        self._images = None
        self._nImages = None

    @property
    def nImages(self):
        """get nImage"""
        return self._nImages

    @property
    def nl(self):
        """Get Neighbor Lists"""
        return self._nl

    @nl.setter
    def nl(self, nl):
        """Set Neighbor Lists"""
        if isinstance(nl, list):
            assert len(nl) == self.nImages
            self._nl = nl
        else:
            self._nl = [ nl ]

        self._delete_derived_info()

    @nl.deleter
    def nl(self):
        """Delete Neighbor Lists"""
        self._nl = None

    def _get_all_x(self, distance):
        """Helper function to get bonds, angles, dihedrals"""
        maxIter = self.nImages
        if len(self.nl) == 1:
            maxIter = 1

        xList = []
        for i in range(maxIter):
            xList.append(get_distance_indices(self.distance_matrix[i], distance))

        return xList

    @property
    def all_bonds(self):
        """Get all Bonds

        Returns a list with indices of bonded atoms for each neighborlist in *self*.
        Atom i is connected to all atoms inside result[i]. Duplicates from PBC are
        removed.
        """
        if self._allBonds is None:
            self._allBonds = self._get_all_x(1)

        return self._allBonds

    @all_bonds.deleter
    def all_bonds(self):
        self._allBonds = None


    @property
    def all_angles(self):
        """Get all angles

        Returns a list with indices of atoms in angles for each neighborlist in *self*.
        Atom i forms an angle to the atoms inside the tuples in result[i]:
        i -- result[i][x][0] -- result[i][x][1]
        where x is in range(number of angles from i).
        """
        if self._allAngles is None:
            self._allAngles = []
            distList = self._get_all_x(2)

            for imI in range(len(distList)):
                self._allAngles.append([])
                #iterate over second neighbors of all atoms
                for iAtom, secNeighs in enumerate(distList[imI]):
                    self._allAngles[-1].append([])
                    if len(secNeighs) == 0:
                        continue
                    firstNeighs = self.all_bonds[imI][iAtom]
                    #iterate over second neighbors of iAtom
                    for kAtom in secNeighs:
                        relevantFirstNeighs = [ idx for idx in firstNeighs if kAtom in self.all_bonds[imI][idx] ]
                        #iterate over all atoms that are connected to iAtom and kAtom
                        for jAtom in relevantFirstNeighs:
                            self._allAngles[-1][-1].append((jAtom, kAtom))

        return self._allAngles

    @all_angles.deleter
    def all_angles(self):
        self._allAngles = None


    @property
    def all_dihedrals(self):
        """Get all dihedrals

        Returns a list with indices of atoms in dihedrals for each neighborlist in *self*.
        Atom i forms a dihedral to the atoms inside the tuples in result[i]:
        i -- result[i][x][0] -- result[i][x][1] -- result[i][x][2]
        where x is in range(number of dihedrals from i).
        """
        if self._allDihedrals is None:
            self._allDihedrals = []
            distList = self._get_all_x(3)

            for imI in range(len(distList)):
                self._allDihedrals.append([])
                for iAtom, thirdNeighs in enumerate(distList[imI]):
                    self._allDihedrals[-1].append([])
                    if len(thirdNeighs) == 0:
                        continue
                    anglesI = self.all_angles[imI][iAtom]
                    #iterate over third neighbors of iAtom
                    for lAtom in thirdNeighs:
                        secondNeighs = [ angle[-1] for angle in anglesI ]
                        firstNeighs = [ angle[0] for angle in anglesI ]
                        relevantSecondNeighs = [ idx for idx in secondNeighs if lAtom in self.all_bonds[imI][idx] ]
                        relevantFirstNeighs = [ firstNeighs[secondNeighs.index(idx)] for idx in relevantSecondNeighs ]
                        #iterate over all atoms that are connected to iAtom and lAtom
                        for jAtom, kAtom in zip(relevantFirstNeighs, relevantSecondNeighs):
                            #remove dihedrals in circles
                            tupl = (jAtom, kAtom, lAtom)
                            if len(set((iAtom, ) + tupl)) != 4:
                                continue
                            #avoid duplicates
                            elif tupl in self._allDihedrals[-1][-1]:
                                continue
                            elif iAtom in tupl:
                                raise RuntimeError("Something is wrong in analysis.all_dihedrals!")
                            self._allDihedrals[-1][-1].append((jAtom, kAtom, lAtom))

        return self._allDihedrals

    @all_dihedrals.deleter
    def all_dihedrals(self):
        self._allDihedrals = None

    @property
    def adjacency_matrix(self):
        """Get the adjacency matrix.

        If not already done, build a list of adjacency matrices for all `self.nl`.
        """

        if self._adjacencyMatrix is None:
            self._adjacencyMatrix = []
            for i in range(len(self.nl)):
                self._adjacencyMatrix.append(self.nl[i].get_connectivity_matrix())

        return self._adjacencyMatrix

    @adjacency_matrix.deleter
    def adjacency_matrix(self):
        self._adjacencyMatrix = None

    @property
    def distance_matrix(self):
        """Get the distance matrix.

        If not already done, build a list of distance matrices for all `self.nl`.
        """

        if self._distanceMatrix is None:
            self._distanceMatrix = []
            for i in range(len(self.nl)):
                self._distanceMatrix.append(get_distance_matrix(self.adjacency_matrix[i]))

        return self._distanceMatrix

    @distance_matrix.deleter
    def distance_matrix(self):
        self._distanceMatrix = None


    @property
    def unique_bonds(self):
        """Get unique bonds.

        Get *all_bonds* i-j without j-i. This is the upper triangle of the
        connectivity matrix (i,j), `i < j`
        """
        bonds = []
        for imI in range(len(self.all_bonds)):
            bonds.append([])
            for iAtom, bonded in enumerate(self.all_bonds[imI]):
                bonds[-1].append([ jAtom for jAtom in bonded if jAtom > iAtom ])

        return bonds


    def _filter_unique(self, l):
        """Helper function to filter for unique lists in a list
        that also contains the reversed items.
        """
        r = []
        #iterate over images
        for imI in range(len(l)):
            r.append([])
            #iterate over atoms
            for i, tuples in enumerate(l[imI]):
                #add the ones where i is smaller than the last element
                r[-1].append([ x for x in tuples if i < x[-1]  ])
        return r

    @property
    def unique_angles(self):
        """Get unique angles.

        Get *all_angles* i-j-k without k-j-i.
        """
        return self._filter_unique(self.all_angles)

    @property
    def unique_dihedrals(self):
        """Get unique dihedrals.

        Get *all_dihedrals* i-j-k-l without l-k-j-i.
        """
        return self._filter_unique(self.all_dihedrals)

    def _get_symbol_idxs(self, imI, sym):
        """Get list of indices of element *sym*"""
        if isinstance(imI, int):
            return [ idx for idx in range(len(self.images[imI])) if self.images[imI][idx].symbol == sym  ]
        else:
            return [ idx for idx in range(len(imI)) if imI[idx].symbol == sym ]

    def _idxTuple2SymbolTuple(self, imI, tup):
        """Converts a tuple of indices to their symbols"""
        return ( self.images[imI][idx].symbol for idx in tup )

    def get_bonds(self, A, B, unique=True):
        """Get bonds from element A to element B"""
        r = []
        for imI in range(len(self.all_bonds)):
            r.append([])
            aIdxs = self._get_symbol_idxs(imI, A)
            if A != B:
                bIdxs = self._get_symbol_idxs(imI, B)
            for idx in aIdxs:
                bonded = self.all_bonds[imI][idx]
                if A == B:
                    r[-1].extend([ (idx, x) for x in bonded if ( x in aIdxs ) and ( x > idx ) ])
                else:
                    r[-1].extend([ (idx, x) for x in bonded if x in bIdxs ])

            if not unique:
                r[-1] +=  [ x[::-1] for x in r[-1] ]

        return r


    def get_angles(self, A, B, C, unique=True):
        """Get angles from given elements A-B-C.

        *B* will be the central atom. The order of the returned indices
        will be A-B-C, if `unique=False` A-B-C AND C-B-A both be returned.
        """
        from itertools import product, combinations, permutations
        r = []
        for imI in range(len(self.all_angles)):
            r.append([])
            #Middle Atom is fixed
            bIdxs = self._get_symbol_idxs(imI, B)
            for bIdx in bIdxs:
                bondedA = [ idx for idx in self.all_bonds[imI][bIdx] if self.images[imI][idx].symbol == A ]
                if len(bondedA) == 0:
                    continue

                if A != C:
                    bondedC = [ idx for idx in self.all_bonds[imI][bIdx] if self.images[imI][idx].symbol == C ]
                    if len(bondedC) == 0:
                        continue

                if A == C:
                    extend = [ (x[0], bIdx, x[1]) for x in list(combinations(bondedA, 2)) ]
                else:
                    extend = list(product(bondedA, [bIdx], bondedC))

                if not unique:
                    extend += [ x[::-1] for x in extend ]

                r[-1].extend(extend)
        return r


    def get_dihedrals(self, A, B, C, D, unique=True):
        """Get dihedrals A-B-C-D.

        If `unique=False` A-B-C-D and D-C-B-A will be returned.
        """
        r = []
        for imI in range(len(self.all_dihedrals)):
            r.append([])
            #get indices of elements
            aIdxs = self._get_symbol_idxs(imI, A)
            bIdxs = self._get_symbol_idxs(imI, B)
            cIdxs = self._get_symbol_idxs(imI, C)
            dIdxs = self._get_symbol_idxs(imI, D)
            for aIdx in aIdxs:
                dihedrals = [ (aIdx, ) + d for d in self.all_dihedrals[imI][aIdx] if ( d[0] in bIdxs ) and ( d[1] in cIdxs ) and ( d[2] in dIdxs ) ]
                if not unique:
                    dihedrals += [ d[::-1] for d in dihedrals ]
                r[-1].extend(dihedrals)

        return r


    def get_bond_value(self, imIdx, idxs, **kwargs):
        """Get bond length idxs[0]-idxs[1] from image imIdx.

        *kwargs* are passed on to :func:`ase.Atoms.get_distance`
        """
        return self.images[imIdx].get_distance(idxs[0], idxs[1], **kwargs)

    def get_angle_value(self, imIdx, idxs, **kwargs):
        """Get angle idxs[0]-idxs[1]-idxs[2] from image imIdx.

        *kwargs* are passed on to :func:`ase.Atoms.get_angle`
        """
        return self.images[imIdx].get_angle(idxs[0], idxs[1], idxs[2], **kwargs)

    def get_dihedral_value(self, imIdx, idxs, **kwargs):
        """Get dihedral idxs[0]-idxs[1]-idxs[2]-idxs[3] from image imIdx.

        *kwargs* are passed on to :func:`ase.Atoms.get_dihedral`
        """
        return self.images[imIdx].get_dihedral(idxs[0], idxs[1], idxs[2], idxs[3], **kwargs)

    def get_values(self, inputList, imageIdx=None, **kwargs):
        """Get Bond/Angle/Dihedral values.

        *inputList* can be any list provided by :meth:`~ase.geometry.analysis.Analysis.get_bonds`,
        :meth:`~ase.geometry.analysis.Analysis.get_angles` or
        :meth:`~ase.geometry.analysis.Analysis.get_dihedrals`.

        Using *imageIdx* (can be integer or slice) the analyzed frames can be specified.
        If *imageIdx* is None, all frames will be analyzed. See :func:`~ase.geometry.analysis.Analysis._get_slice` for details.

        *kwargs* is passed on to the :class:`~ase.Atoms` classes functions for
        retrieving the values.

        The type of value requested is determined from the length of the tuple inputList[0][0].
        The methods from the :class:`~ase.Atoms` class are used.
        """

        sl = self._get_slice(imageIdx)

        #get method to call from length of inputList
        if len(inputList[0][0]) == 2:
            get = self.get_bond_value
        elif len(inputList[0][0]) == 3:
            get = self.get_angle_value
        elif len(inputList[0][0]) == 4:
            get = self.get_dihedral_value
        else:
            raise ValueError("inputList in ase.geometry.analysis.Analysis.get_values has a bad shape.")

        #check if length of slice and inputList match
        singleNL = False
        if len(inputList) != len(self.images[sl]):
            #only one nl for all images
            if len(inputList) == 1 and len(self.nl) == 1:
                singleNL = True
            else:
                raise RuntimeError("Length of inputList does not match length of \
                        images requested, but it also is not one item long.")

        r = []
        for inputIdx, image in enumerate(self.images[sl]):
            imageIdx = self.images.index(image)
            r.append([])
            #always use first list from input if only a single neighborlist was used
            if singleNL:
                inputIdx = 0
            for tupl in inputList[inputIdx]:
                r[-1].append(get(imageIdx, tupl, **kwargs))

        return r

    def get_rdf(self, rmax, nbins, imageIdx=None, elements=None, return_dists=False):
        """Get RDF.

        Wrapper for :meth:`ase.ga.utilities.get_rdf` with more selection possibilities.

        * `rmax` (type: float): Maximum distance of RDF.
        * `nbins` (type: int): Number of bins to devide RDF.
        * `imageIdx` (type: int/slice/None): Images to analyze, see :func:`~ase.geometry.analysis.Analysis._get_slice` for details.
        * `elements` (type: str/int/list/tuple): Make partial RDFs.

        If elements is `None`, a full RDF is calculated. If elements is an *integer* or a *list/tuple
        of integers*, only those atoms will contribute to the RDF (like a mask). If elements
        is a *string* or a *list/tuple of strings*, only Atoms of those elements will contribute.

        """

        sl = self._get_slice(imageIdx)

        r = []
        el = None

        for image in self.images[sl]:
            if elements is None:
                tmpImage = image
            #integers
            elif isinstance(elements, int):
                tmpImage = Atoms(cell=image.get_cell(), pbc=image.get_pbc())
                tmpImage.append(image[elements])
            #strings
            elif isinstance(elements, str):
                tmpImage = Atoms(cell=image.get_cell(), pbc=image.get_pbc())
                for idx in self._get_symbol_idxs(image, elements):
                    tmpImage.append(image[idx])
            #lists
            elif isinstance(elements, list) or isinstace(elements, tuple):
                #list of ints
                if all(isinstance(x, int) for x in elements):
                    if len(elements) == 2:
                        #use builtin get_rdf mask
                        el = elements
                        tmpImage = image
                    else:
                        #create dummy image
                        tmpImage = Atoms(cell=image.get_cell(), pbc=image.get_pbc())
                        for idx in elements:
                            tmpImage.append(image[idx])
                #list of strings
                elif all(isinstance(x, str) for x in elements):
                    tmpImage = Atoms(cell=image.get_cell(), pbc=image.get_pbc())
                    for element in elements:
                        for idx in self._get_symbol_idxs(image, element):
                            tmpImage.append(image[idx])
                else:
                    raise ValueError("Unsupported type of elements given in ase.geometry.analysis.Analysis.get_rdf!")
            else:
                raise ValueError("Unsupported type of elements given in ase.geometry.analysis.Analysis.get_rdf!")

            r.append(get_rdf(tmpImage, rmax, nbins, elements=el, no_dists=(not return_dists)))
        return r



