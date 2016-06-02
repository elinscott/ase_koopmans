from ase.atoms import Atoms
from ase.quaternions import Quaternions
from ase.calculators.singlepoint import SinglePointCalculator
from ase.parallel import paropen
from ase.utils import basestring
from collections import deque
import numpy as np


def read_lammps_dump(fileobj, index=-1, order=True, atomsobj=Atoms):
    """Method which reads a LAMMPS dump file.

    order: Order the particles according to their id. Might be faster to
    switch it off.
    """
    if isinstance(fileobj, basestring):
        f = paropen(fileobj)
    else:
        f = fileobj

    # load everything into memory
    lines = deque(f.readlines())

    n_atoms = 0
    images = []

    while len(lines) > n_atoms:
        line = lines.popleft()

        if 'ITEM: TIMESTEP' in line:

            n_atoms = 0
            lo = []
            hi = []
            tilt = []

        if 'ITEM: NUMBER OF ATOMS' in line:
            line = lines.popleft()
            n_atoms = int(line.split()[0])
            
        if 'ITEM: BOX BOUNDS' in line:
            # save labels behind "ITEM: BOX BOUNDS" in triclinic case
            # (>=lammps-7Jul09)
            tilt_items = line.split()[3:]
            for i in range(3):
                line = lines.popleft()
                fields = line.split()
                lo.append(float(fields[0]))
                hi.append(float(fields[1]))
                if (len(fields) >= 3):
                    tilt.append(float(fields[2]))

            # determine cell tilt (triclinic case!)
            if len(tilt) >= 3:
                # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS"
                # to assign tilt (vector) elements ...
                if len(tilt_items) >= 3:
                    xy = tilt[tilt_items.index('xy')]
                    xz = tilt[tilt_items.index('xz')]
                    yz = tilt[tilt_items.index('yz')]
                # ... otherwise assume default order in 3rd column
                # (if the latter was present)
                else:
                    xy = tilt[0]
                    xz = tilt[1]
                    yz = tilt[2]
            else:
                xy = xz = yz = 0
            xhilo = (hi[0] - lo[0]) - abs(xy) - abs(xz)
            yhilo = (hi[1] - lo[1]) - abs(yz)
            zhilo = (hi[2] - lo[2])
            celldispx = lo[0] - min(0, xy) - min(0, xz)
            celldispy = lo[1] - min(0, yz)
            celldispz = lo[2]

            cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]
            celldisp = [[celldispx, celldispy, celldispz]]
        
        if 'ITEM: ATOMS' in line:
            colnames = line.split()[2:]
            datarows = [lines.popleft() for _ in range(n_atoms)]
            # lines = lines[n_atoms:]
            integer_cols = (colnames.index('id'), colnames.index('type'))
            ids, types = np.loadtxt(datarows, usecols=integer_cols, dtype=int).T
            data = np.loadtxt(datarows)
            if order:
                data = data[ids-1, :]
                types = types[ids-1]
            
            def get_quantity(labels):
                try:
                    cols = [colnames.index(label) for label in labels]
                    return data[:, cols]
                except ValueError:
                    return []
            
            positions = get_quantity(['x', 'y', 'z'])
            scaled_positions = get_quantity(['xs', 'ys', 'zs'])
            velocities = get_quantity(['vx', 'vy', 'vz'])
            charges = get_quantity(['q'])[:]
            forces = get_quantity(['fx', 'fy', 'fz'])
            quaternions = get_quantity(['c_q[1]', 'c_q[2]', 'c_q[3]', 'c_q[4]'])
            
            if len(quaternions):
                images.append(Quaternions(symbols=types,
                                          positions=positions,
                                          cell=cell, celldisp=celldisp,
                                          quaternions=quaternions))
            elif len(positions):
                images.append(atomsobj(
                    symbols=types, positions=positions,
                    celldisp=celldisp, cell=cell))
            elif len(scaled_positions):
                images.append(atomsobj(
                    symbols=types, scaled_positions=scaled_positions,
                    celldisp=celldisp, cell=cell))

            if len(velocities):
                images[-1].set_velocities(velocities)
            if len(charges):
                images[-1].set_initial_charges(charges)
            if len(forces):
                calculator = SinglePointCalculator(images[-1],
                                                   energy=0.0, forces=forces)
                images[-1].set_calculator(calculator)

    return images[index]
