import gzip
from os.path import splitext
from collections import deque
import numpy as np

from ase.atoms import Atoms
from ase.quaternions import Quaternions
from ase.calculators.singlepoint import SinglePointCalculator
from ase.parallel import paropen
from ase.utils import basestring


def read_lammps_dump(infileobj, **kwargs):
    """Method which reads a LAMMPS dump file.

       LAMMPS chosen output method depends on the chosen suffix:
        - .bin  : binary file
        - .gz   : output piped through gzip
        - .mpiio: using mpiio (should be like cleartext, with different ordering)
        - else  : normal clear-text format

    :param infileobj: string to file, opened file or file-like stream

    """
    # !TODO: add support for lammps-regex expression naming schemes
    # (output per processor and timestep wildcards)

    if isinstance(infileobj, basestring):
        suffix = splitext(infileobj)[-1]
        if suffix == ".bin":
            fileobj = paropen(infileobj, 'rb')
        if suffix == ".gz":
            # !TODO: save for parallel execution?
            fileobj = gzip.open(infileobj, 'rb')
        else:
            fileobj = paropen(infileobj)
    else:
        suffix = splitext(infileobj.__name__)
        fileobj = infileobj

    if suffix == ".bin":
        return read_lammps_dump_bin(fileobj, **kwargs)

    return read_lammps_dump_string(fileobj, **kwargs)


def lammps_data_to_ase_atoms(data, colnames, cell, celldisp, pbc = False,
                             atomsobj=Atoms, order=True, specorder=None):
    """Extract positions and other per-atom parameters and create Atoms

    :param datarows: per atom data
    :param colnames: index for data
    :param cell: cell dimision
    :param celldisp: origin shift
    :param pbc: periodic boundaries
    :param atomsobj: Class to create Atoms object
    :param order: sort atoms by id
    :param specorder: map lammps types to species
    :returns: Atoms object
    :rtype: Atoms

    """
    ids = data[:, colnames.index('id')].astype(int)
    types = data[:, colnames.index('type')].astype(int)
    if order:
        data = data[ids-1, :]
        types = types[ids-1]

    # reconstruct types from given specorder
    if specorder:
        types = [specorder[t-1] for t in types]

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
        out_atoms = Quaternions(symbols=types,
                                positions=positions,
                                cell=cell, celldisp=celldisp,
                                quaternions=quaternions)
    elif len(positions):
        out_atoms = atomsobj(symbols=types, positions=positions,
                             celldisp=celldisp, cell=cell)
    elif len(scaled_positions):
        out_atoms = atomsobj(symbols=types, scaled_positions=scaled_positions,
                             celldisp=celldisp, cell=cell)

    if len(velocities):
        out_atoms.set_velocities(velocities)
    if len(charges):
        out_atoms.set_initial_charges(charges)
    if len(forces):
        calculator = SinglePointCalculator(out_atoms,
                                           energy=0.0, forces=forces)
        out_atoms.set_calculator(calculator)

    return out_atoms


def read_lammps_dump_string(fileobj, index=-1, order=True,
                            atomsobj=Atoms, specorder=None):
    """Process cleartext lammps dumpfiles

    :param fileobj: filestream providing the trajectory data
    :param index: if containing multiple images, which to return (default: the last)
    :param order: Order the particles according to their id. Might be faster to
    switch it off.
    :param atomsobj: function to create ase-Atoms object 
    :param specorder: list of species in data-file 
    (usually .dump files to not contain type to species mapping)
    :returns: list of Atoms objects
    :rtype: list
    """

    # load everything into memory
    lines = deque(fileobj.readlines())

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
            for _ in range(3):
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
            data = np.loadtxt(datarows)
            out_atoms = lammps_data_to_ase_atoms(
                data=data, colnames=colnames, cell=cell, celldisp=celldisp,
                atomsobj=Atoms, order=order, specorder=specorder)
            images.append(out_atoms)

    return images[index]


def read_lammps_dump_bin(fileobj):
    raise NotImplementedError
