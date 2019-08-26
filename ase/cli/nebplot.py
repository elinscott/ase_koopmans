import sys

import numpy as np

from ase.db import connect
from ase.build import bulk
from ase.io import read, write
from ase.visualize import view
from ase.build import molecule
from ase.atoms import Atoms
from ase.symbols import string2symbols
from ase.data import ground_state_magnetic_moments
from ase.data import atomic_numbers, covalent_radii
from ase.neb import NEBTools


class CLICommand:
    """Analyze NEB trajectories by making band plots.

    One file:

        ase nebplot neb.traj

    Multiple files:

        ase nebplot neb1.traj neb2.traj

    Specify output:

        ase nebplot neb1.traj neb2.traj myfile.pdf
    """

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('filenames', nargs='+', help='one or more trajectory files to analyze')
        add('output', nargs='?', help='optional name of output file, default=nebplots.pdf')
        add('-i', dest='n_images', default=None,
            help='number of images per band, guessed if not supplied')
        add('-x', dest='constant_x', action='store_true',
            help='use a single x axis scale for all plots')
        add('-y', dest='constant_y', action='store_true',
            help='use a single y axis scale for all plots')

    @staticmethod
    def run(args, parser):
        # Nothing will ever be stored in args.output; need to manually find
        # if its supplied by checking extensions.
        if args.filenames[-1].endswith('.pdf'):
            args.output = args.filenames.pop(-1)

        nebtools = NEBTools(images=args.filenames)
        if args.output:
            nebtools.plot_bands(constant_x=args.constant_x,
                                constant_y=args.constant_y,
                                nimages=args.n_images,
                                label=args.output[:-4])
        else:
            nebtools.plot_bands(constant_x=args.constant_x,
                                constant_y=args.constant_y,
                                nimages=args.n_images)



