import numpy as np
from ase.build.general_surface import surface
from ase.atoms import string2symbols


def find_z_layers(atoms, tolerance=0.2, max_width=0.4):
    """
    Find layers, returning a list of arrays of atom indices for each layer.

    For example, atoms[layers[i]] will be the ith layer."""

    args = np.argsort(atoms.positions[:, 2])
    zcoords = atoms.positions[:, 2][args]
    layers = []

    def add_layer(layer):
        layer_width = zcoords[layer[-1]] - zcoords[layer[0]]
        if layer_width > max_width:
            raise ValueError('Layers not well defined: width={} > '
                             'max_width={}'.format(layer_width, max_width))
        layers.append(layer)

    layer = [0]
    for i in range(1, len(atoms)):
        if zcoords[i] - zcoords[i - 1] > tolerance:
            add_layer(layer)
            layer = []
        layer.append(i)
    add_layer(layer)
    return [args[layer] for layer in layers]


def Smart_Surface(lattice, indices, layers, vacuum=None, tol=1e-10,
                  termination=None, return_all=False, verbose=False):
    """
    Create surface from a given lattice and Miller indices. Smart Surface
    allows you to ask for a certain atomic termination.

    lattice: Atoms object or str
        Bulk lattice structure of alloy or pure metal.  Note that the
        unit-cell must be the conventional cell - not the primitive cell.
        One can also give the chemical symbol as a string, in which case the
        correct bulk lattice will be generated automatically.
    indices: sequence of three int
        Surface normal in Miller indices (h,k,l).
    layers: int
        Number of equivalent layers of the slab. (not the same as the layers
        you choose from for terminations)
    vacuum: float
        Amount of vacuum added on both sides of the slab.
    termination: str
        the atoms you wish to be in the top layer. There may be many such
        terminations, this function returns all terminations with the same
        atomic composition.
        e.g. 'O' will return oxygen terminated surfaces.
        e.g.'TiO' will return surfaces terminated with layers containing both O
        and Ti

    """
    lats = translate_lattice(lattice, indices)
    return_surfs = []
    check = []
    check2 = []
    for item in lats:
        check3 = True
        surf = surface(item, indices, layers, vacuum=vacuum, tol=tol)
        positions = surf.get_scaled_positions().flatten()
        for i, value in enumerate(positions):
            if value >= 1 - tol:  # move things closer to zero within tol
                positions[i] -= 1
        surf.set_scaled_positions(np.reshape(positions, (len(surf), 3)))
        rep = find_z_layers(surf)
        if termination is not None:
            comp = [surf.get_chemical_symbols()[a] for a in rep[-1]]
            term = string2symbols(termination)
            # list atoms in top layer and not in requested termination
            check = [a for a in comp if a not in term]
            # list of atoms in requested termination and not in top layer
            check2 = [a for a in term if a not in comp]
        if len(return_surfs) > 0:
            pos_diff = [a.get_positions() - surf.get_positions()
                        for a in return_surfs]
            for i, su in enumerate(pos_diff):
                similarity_test = su.flatten() < tol * 1000
                if similarity_test.all():
                    # checks if surface is too similar to another surface
                    check3 = False
        if not check3:
            continue
        if return_all is True:
            pass
        elif check != [] or check2 != [] and check3:
            continue
        return_surfs.append(surf)
    return return_surfs

            
def translate_lattice(lattice, indices, tol=10**-3):
    """
    translates a bulk unit cell along a normal vector given by the a set of
    miller indices to the next symetric position. This is used to control the
    termination of the surface in the smart_surface command
    inputs:
        lattice: Atoms object
            atoms object of the bulk unit cell
        indices: 1x3 list,tuple, or numpy array
            the miller indices you wish to cut along.
    returns:
        lattice_list: list of Atoms objects
            a list of all the different translations of the unit cell that will
            yield different terminations of a surface cut along the miller
            indices provided.
    """
    lattice_list = []
    cell = lattice.get_cell()
    pt = [0, 0, 0]
    h, k, l = indices
    millers = list(indices)
    for index, item in enumerate(millers):
        if item == 0:
            millers[index] = 10**9  # make zeros large numbers
        elif pt == [0, 0, 0]:  # for numerical stability
            pt = list(cell[index] / float(item) / np.linalg.norm(cell[index]))
    h1, k1, l1 = millers
    N = np.array(cell[0] / h1 + cell[1] / k1 + cell[2] / l1)
    n = N / np.linalg.norm(N)  # making a unit vector normal to cut plane
    # finding distance from cut plan vector
    d = [np.round(np.dot(n, (a - pt)) * n, 5) for
         a in lattice.get_scaled_positions()]
    duplicates = []
    for i, item in enumerate(d):
        g = [True for a in d[i + 1:] if np.linalg.norm(a - item) < tol]
        if g != []:
            duplicates.append(i)
    duplicates.reverse()
    for i in duplicates:
        del d[i]
    # put distance to the plane at the end of the array
    for i, item in enumerate(d):
        d[i] = np.append(item,
                         np.dot(n, (lattice.get_scaled_positions()[i] - pt)))
    d = np.array(d)
    d = d[d[:, 3].argsort()]  # sort by distance to the plane
    d = [a[:3] for a in d]  # remove distance
    d = list(d)  # make it a list again
    for i in d:
        """
        The above method gives you the boundries of between terminations that
        will allow you to build a complete set of terminations. However, it
        does not return all the boundries. Thus you must check both above and
        below the boundry, and not stray too far from the boundry. If you move
        too far away, you risk hitting another boundry you did not find.
        """
        lattice1 = lattice.copy()
        displacement = (h * cell[0] + k * cell[1] + l * cell[2]) \
            * (i + 10 ** -8)
        lattice1.positions -= displacement
        lattice_list.append(lattice1)
        lattice1 = lattice.copy()
        displacement = (h * cell[0] + k * cell[1] + l * cell[2]) \
            * (i - 10 ** -8)
        lattice1.positions -= displacement
        lattice_list.append(lattice1)
    return lattice_list
