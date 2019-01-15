import numpy as np
from ase.geometry import find_mic


def fit_raw(energies, forces, positions, cell=None, pbc=None):
    E = energies
    F = forces
    R = positions
    E = np.array(E) - E[0]
    n = len(E)
    Efit = np.empty((n - 1) * 20 + 1)
    Sfit = np.empty((n - 1) * 20 + 1)

    s = [0]
    dR = np.zeros_like(R)
    for i in range(n):
        if i < n - 1:
            dR[i] = R[i + 1] - R[i]
            if cell is not None and pbc is not None:
                dR[i], _ = find_mic(dR[i], cell, pbc)
            s.append(s[i] + np.sqrt((dR[i]**2).sum()))
        else:
            dR[i] = R[i] - R[i - 1]
            if cell is not None and pbc is not None:
                dR[i], _ = find_mic(dR[i], cell, pbc)

    lines = []
    dEds0 = None
    for i in range(n):
        d = dR[i]
        if i == 0:
            ds = 0.5 * s[1]
        elif i == n - 1:
            ds = 0.5 * (s[-1] - s[-2])
        else:
            ds = 0.25 * (s[i + 1] - s[i - 1])

        d = d / np.sqrt((d**2).sum())
        dEds = -(F[i] * d).sum()
        x = np.linspace(s[i] - ds, s[i] + ds, 3)
        y = E[i] + dEds * (x - s[i])
        lines.append((x, y))

        if i > 0:
            s0 = s[i - 1]
            s1 = s[i]
            x = np.linspace(s0, s1, 20, endpoint=False)
            c = np.linalg.solve(np.array([(1, s0, s0**2, s0**3),
                                          (1, s1, s1**2, s1**3),
                                          (0, 1, 2 * s0, 3 * s0**2),
                                          (0, 1, 2 * s1, 3 * s1**2)]),
                                np.array([E[i - 1], E[i], dEds0, dEds]))
            y = c[0] + x * (c[1] + x * (c[2] + x * c[3]))
            Sfit[(i - 1) * 20:i * 20] = x
            Efit[(i - 1) * 20:i * 20] = y

        dEds0 = dEds

    Sfit[-1] = s[-1]
    Efit[-1] = E[-1]
    return ForceFit(s, E, Sfit, Efit, lines)

from collections import namedtuple

class ForceFit(namedtuple('ForceFit', ['s', 'E', 'Sfit', 'Efit', 'lines'])):
    def plot(self, ax=None):
        import matplotlib.pyplot as plt
        if ax is None:
            ax = plt.gca()

        ax.plot(self.s, self.E, 'o')
        for x, y in self.lines:
            ax.plot(x, y, '-g')
        ax.plot(self.Sfit, self.Efit, 'k-')
        ax.set_xlabel(r'path [$\AA$]')
        ax.set_ylabel('energy [eV]')
        Ef = max(self.E) - self.E[0]
        Er = max(self.E) - self.E[-1]
        dE = self.E[-1] - self.E[0]
        ax.set_title('$E_\mathrm{f} \\approx$ %.3f eV; '
                     '$E_\mathrm{r} \\approx$ %.3f eV; '
                     '$\\Delta E$ = %.3f eV'
                     % (Ef, Er, dE))
        return ax

def fit_images(images):
    R = [atoms.positions for atoms in images]
    E = [atoms.get_potential_energy() for atoms in images]
    F = [atoms.get_forces() for atoms in images]  # XXX force consistent???
    A = images[0].cell
    pbc = images[0].pbc
    return fit_raw(E, F, R, A, pbc)


def energy_force_curve(images, ax=None):
    if ax is None:
        import matplotlib.pyplot as plt
        ax = plt.gca()

    nim = len(images)

    accumulated_distance = 0.0

    energies = [atoms.get_potential_energy()
                for atoms in images]

    for i in range(nim):
        atoms = images[i]
        f_ac = atoms.get_forces()

        if i < nim - 1:
            rightpos = images[i + 1].positions
        else:
            rightpos = atoms.positions

        if i > 0:
            leftpos = images[i - 1].positions
        else:
            leftpos = atoms.positions

        disp_ac = rightpos - leftpos

        def total_displacement(disp):
            disp_a = (disp**2).sum(axis=1)**.5
            return sum(disp_a)

        disp = 0.5 * total_displacement(disp_ac)
        dE_fdotr = -0.5 * np.vdot(f_ac.ravel(), disp_ac.ravel())

        linescale = 0.4
        x1 = accumulated_distance - disp * linescale
        x2 = accumulated_distance + disp * linescale
        y1 = energies[i] - dE_fdotr * linescale
        y2 = energies[i] + dE_fdotr * linescale

        ax.plot([x1, x2], [y1, y2], 'b-')
        ax.plot(accumulated_distance, energies[i], 'bo')
        ax.set_ylabel('Energy [eV]')
        ax.set_xlabel('Accumulative distance [Å]')
        accumulated_distance += total_displacement(rightpos - atoms.positions)
    return ax

if __name__ == '__main__':
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.md import VelocityVerlet
    from ase.units import fs
    from ase.io import read

    atoms = bulk('Au', cubic=True) * (2, 1, 1)
    atoms.calc = EMT()
    atoms.rattle(stdev=0.05)

    md = VelocityVerlet(atoms, timestep=6.0 * fs, trajectory='tmp.traj')
    md.run(steps=50)
    images = read('tmp.traj', ':')

    import matplotlib.pyplot as plt
    ax = energy_force_curve(images)
    plt.savefig('tmp.png')
