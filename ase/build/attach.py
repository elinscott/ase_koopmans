import numpy as np


def nearest(atoms1, atoms2):
    """Return indices of nearest atoms"""
    vd_aac = np.empty((len(atoms1), len(atoms2), 3))
    pos2_ac = atoms2.get_positions()
    for ia, atom in enumerate(atoms1):
        vd_aac[ia] = atom.position - pos2_ac
    d2_aa = (vd_aac**2).sum(axis=2)

    return np.argwhere(d2_aa == d2_aa.min())[0]


def attach(atoms1, atoms2, distance, direction=(1, 0, 0),
           maxiter=100, accuracy=1e-5):
    """Attach two structures

    Parameters:
    atoms1: Atoms
    atoms2: Atoms
    distance: float
      minimal distance (Angstrom)
    direction: unit vector (3 floats)
      relative direction between center of masses
    maxiter: int
      maximal number of iterations to get required distance, default 100
    accuracy: float
      required accuracy for minimal distance (Angstrom), default 1e-5
    """
    direction = np.array(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    assert len(direction) == 3
    dist2 = distance**2
    acc2 = accuracy**2
    
    cm1 = atoms1.get_center_of_mass()
    d1max = np.dot(atoms1.get_positions() - cm1, direction).max()
    cm2 = atoms2.get_center_of_mass()
    d2max = np.dot(cm2 - atoms2.get_positions(), direction).max()

    # first guess
    atoms2.translate(cm1 - cm2 +
                     direction * (distance + d1max + d2max))
    i1, i2 = nearest(atoms1, atoms2)

    for i in range(maxiter):   
        dv_c = atoms2[i2].position - atoms1[i1].position
        dv2 = (dv_c**2).sum()

        # we need to move
        vcost = np.dot(dv_c, direction)
        a = np.sqrt(max(0, dist2 - dv2 + vcost**2))
        atoms2.translate(direction * (a - vcost))

        dv_c = atoms2[i2].position - atoms1[i1].position
        dv2 = (dv_c**2).sum()
        
        if (np.sqrt(dv2) - distance)**2 <= acc2:
            return atoms1 + atoms2

    raise
    
