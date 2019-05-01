from math import pi, sin, cos
import numpy as np


def bz_vertices(icell):
    from scipy.spatial import Voronoi
    I = (np.indices((3, 3, 3)) - 1).reshape((3, 27))
    G = np.dot(icell.T, I).T
    vor = Voronoi(G)
    bz1 = []
    for vertices, points in zip(vor.ridge_vertices, vor.ridge_points):
        if -1 not in vertices and 13 in points:
            normal = G[points].sum(0)
            normal /= (normal**2).sum()**0.5
            bz1.append((vor.vertices[vertices], normal))
    return bz1


def bz_plot(cell, vectors=False, paths=None, points=None,
            elev=None, scale=1, interactive=False,
            pointstyle=None, fig=None, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        if fig is None:
            fig = plt.gcf()
    dimensions = np.sum(cell.pbc)
    assert dimensions > 0, 'No BZ for 0D!'

    if dimensions == 3:
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import proj3d
        from matplotlib.patches import FancyArrowPatch
        Axes3D  # silence pyflakes

        class Arrow3D(FancyArrowPatch):
            def __init__(self, xs, ys, zs, *args, **kwargs):
                FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
                self._verts3d = xs, ys, zs

            def draw(self, renderer):
                xs3d, ys3d, zs3d = self._verts3d
                xs, ys, zs = proj3d.proj_transform(xs3d, ys3d,
                                                   zs3d, renderer.M)
                self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
                FancyArrowPatch.draw(self, renderer)
        azim = pi / 5
        elev = elev or pi / 6
        x = sin(azim)
        y = cos(azim)
        view = [x * cos(elev), y * cos(elev), sin(elev)]
        if ax is None:
            ax = fig.gca(projection='3d')
    elif dimensions == 2:
        # 2d in xy
        assert all(abs(cell[2][0:2]) < 1e-6) and all(abs(cell.T[2]
                                                         [0:2]) < 1e-6)
        # We can add a third dimension for practical purposes which makes
        # the code work much like the 3D case. Then we can just do a normal
        # Voronoi diagram.
        cell = cell.copy()
        cell[2, 2] += 10000.
        ax = plt.axes()
    else:
        # 1d in x
        assert (all(abs(cell[2][0:2]) < 1e-6) and
                all(abs(cell.T[2][0:2]) < 1e-6) and
                abs(cell[0][1]) < 1e-6 and abs(cell[1][0]) < 1e-6)
        ax = plt.axes()
        # Similar to 2D we can add two extra dimensions
        cell = cell.copy()
        cell[1, 1] += 10000.
        cell[2, 2] += 10000.
        ax = plt.axes()

    icell = np.linalg.inv(cell).T
    kpoints = points
    bz1 = bz_vertices(icell)

    maxp = 0.0
    if dimensions == 1:
        x = np.array([-0.5 * icell[0, 0],
                      0.5 * icell[0, 0],
                      -0.5 * icell[0, 0]])
        y = np.array([0, 0, 0])
        ax.plot(x, y, c='k', ls='-')
        maxp = icell[0, 0]
    else:
        for points, normal in bz1:
            x, y, z = np.concatenate([points, points[:1]]).T
            if dimensions == 3:
                if np.dot(normal, view) < 0 and not interactive:
                    ls = ':'
                else:
                    ls = '-'
                    ax.plot(x, y, z, c='k', ls=ls)
            elif dimensions == 2:
                ax.plot(x, y, c='k', ls='-')
            maxp = max(maxp, points.max())

    if vectors:
        if dimensions == 3:
            ax.add_artist(Arrow3D([0, icell[0, 0]],
                                  [0, icell[0, 1]],
                                  [0, icell[0, 2]],
                                  mutation_scale=20, lw=1,
                                  arrowstyle='-|>', color='k'))
            ax.add_artist(Arrow3D([0, icell[1, 0]],
                                  [0, icell[1, 1]],
                                  [0, icell[1, 2]],
                                  mutation_scale=20, lw=1,
                                  arrowstyle='-|>', color='k'))
            ax.add_artist(Arrow3D([0, icell[2, 0]],
                                  [0, icell[2, 1]],
                                  [0, icell[2, 2]],
                                  mutation_scale=20, lw=1,
                                  arrowstyle='-|>', color='k'))
            maxp = max(maxp, 0.6 * icell.max())
        elif dimensions == 2:
            ax.arrow(0, 0, icell[0, 0], icell[0, 1],
                     lw=1, color='k',
                     length_includes_head=True,
                     head_width=0.03, head_length=0.05)
            ax.arrow(0, 0, icell[1, 0], icell[1, 1],
                     lw=1, color='k',
                     length_includes_head=True,
                     head_width=0.03, head_length=0.05)
            maxp = max(maxp, icell.max())
        else:
            ax.arrow(0, 0, icell[0, 0], 0,
                     lw=1, color='k',
                     length_includes_head=True,
                     head_width=0.03, head_length=0.05)
            maxp = max(maxp, icell.max())

    if paths is not None:
        for names, points in paths:
            x, y, z = np.array(points).T
            if dimensions == 3:
                ax.plot(x, y, z, c='r', ls='-', marker='.')
            elif dimensions in [1, 2]:
                ax.plot(x, y, c='r', ls='-')
                
            for name, point in zip(names, points):
                x, y, z = point
                if name == 'G':
                    name = '\\Gamma'
                elif len(name) > 1:
                    import re
                    m = re.match(r'^(\D+?)(\d*)$', name)
                    if m is None:
                        raise ValueError('Bad label: {}'.format(name))
                    name, num = m.group(1, 2)
                    if num:
                        name = '{}_{{{}}}'.format(name, num)
                if dimensions == 3:
                    ax.text(x, y, z, '$' + name + '$',
                            ha='center', va='bottom', color='r')
                elif dimensions == 2:
                    if abs(z) < 1e-6:
                        ax.text(x, y, '$' + name + '$',
                                ha='center', va='bottom', color='r')
                else:
                    if abs(y) < 1e-6 and abs(z) < 1e-6:
                        ax.text(x, y, '$' + name + '$',
                                ha='center', va='bottom', color='r')

    if kpoints is not None:
        kw = {'c': 'b'}
        if pointstyle is not None:
            kw.update(pointstyle)
        for p in kpoints:
            if dimensions == 3:
                ax.scatter(p[0], p[1], p[2], **kw)
            elif dimensions == 2:
                ax.scatter(p[0], p[1], c='b')
            else:
                ax.scatter(p[0], 0, c='b')

    ax.set_axis_off()
    ax.autoscale_view(tight=True)
    if dimensions == 3:
        s = maxp / 0.5 * 0.45 * scale
    elif dimensions in [1, 2]:
        s = maxp * 1.05
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)
    if dimensions == 3:
        ax.set_zlim(-s, s)
    ax.set_aspect('equal')
    if dimensions == 3:
        ax.view_init(azim=azim / pi * 180, elev=elev / pi * 180)
    return ax
