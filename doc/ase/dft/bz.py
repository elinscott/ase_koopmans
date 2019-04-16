# creates: bztable.rst
import numpy as np
import matplotlib.pyplot as plt

from ase.dft.kpoints import (get_special_points, special_paths,
                             parse_path_string)
from ase.dft.bz import bz3d_plot

from ase.geometry.bravais import all_variants


header = """\

Brillouin zone data
-------------------

.. list-table::
    :widths: 10 15 45
"""


entry = """\
    * - {name} ({longname})
      - {bandpath}
      - .. image:: {fname}
            :width: 40 %
"""

with open('bztable.rst', 'w') as fd:
    print(header, file=fd)

    for i, lat in enumerate(all_variants()):
        id = '{:02d}.{}'.format(i, lat.variant.name)
        imagefname = '{}.svg'.format(id)
        txt = entry.format(name=lat.variant.name,
                           longname=lat.longname,
                           bandpath=lat.bandpath().labelseq,
                           fname=imagefname)
        print(txt, file=fd)
        ax = lat.plot_bz()
        fig = ax.get_figure()
        fig.savefig(imagefname, bbox_inches='tight')
        fig.clear()
