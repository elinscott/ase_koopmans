import os
import subprocess
import sys
import tempfile

from ase_koopmans.io import write
import ase_koopmans.parallel as parallel


def _pipe_to_gui(atoms, repeat, block):
    from io import BytesIO
    buf = BytesIO()
    write(buf, atoms, format='traj')

    args = [sys.executable, '-m', 'ase_koopmans', 'gui', '-']
    if repeat:
        args.append(' --repeat={},{},{}'.format(*repeat))

    proc = subprocess.Popen(args,
                            stdin=subprocess.PIPE)
    proc.stdin.write(buf.getvalue())
    proc.stdin.close()
    if block:
        proc.wait()


def view(atoms, data=None, viewer='ase_koopmans', repeat=None, block=False):
    # Ignore for parallel calculations:
    if parallel.world.size != 1:
        return

    vwr = viewer.lower()

    if vwr == 'ase_koopmans':
        _pipe_to_gui(atoms, repeat, block)
        return

    if vwr == 'vmd':
        format = 'cube'
        command = 'vmd'
    elif vwr == 'rasmol':
        format = 'proteindatabank'
        command = 'rasmol -pdb'
    elif vwr == 'xmakemol':
        format = 'xyz'
        command = 'xmakemol -f'
    elif vwr == 'gopenmol':
        format = 'xyz'
        command = 'rungOpenMol'
    elif vwr == 'avogadro':
        format = 'cube'
        command = 'avogadro'
    elif vwr == 'sage':
        from ase_koopmans.visualize.sage import view_sage_jmol
        view_sage_jmol(atoms)
        return
    elif vwr in ('ngl', 'nglview'):
        from ase_koopmans.visualize.ngl import view_ngl
        return view_ngl(atoms)
    elif vwr == 'x3d':
        from ase_koopmans.visualize.x3d import view_x3d
        return view_x3d(atoms)
    elif vwr == 'paraview':
        # macro for showing atoms in paraview
        macro = """\
from paraview.simple import *
version_major = servermanager.vtkSMProxyManager.GetVersionMajor()
source = GetActiveSource()
renderView1 = GetRenderView()
atoms = Glyph(Input=source,
              GlyphType='Sphere',
#              GlyphMode='All Points',
              Scalars='radii',
              ScaleMode='scalar',
              )
RenameSource('Atoms', atoms)
atomsDisplay = Show(atoms, renderView1)
if version_major <= 4:
    atoms.SetScaleFactor = 0.8
    atomicnumbers_PVLookupTable = GetLookupTableForArray( "atomic numbers", 1)
    atomsDisplay.ColorArrayName = ('POINT_DATA', 'atomic numbers')
    atomsDisplay.LookupTable = atomicnumbers_PVLookupTable
else:
    atoms.ScaleFactor = 0.8
    ColorBy(atomsDisplay, 'atomic numbers')
    atomsDisplay.SetScalarBarVisibility(renderView1, True)
Render()
        """
        script_name = os.path.join(tempfile.gettempdir(), 'draw_atoms.py')
        with open(script_name, 'w') as f:
            f.write(macro)
        format = 'vtu'
        command = 'paraview --script=' + script_name
    else:
        raise RuntimeError('Unknown viewer: ' + viewer)

    fd, filename = tempfile.mkstemp('.' + format, 'ase_koopmans-')
    if repeat is not None:
        atoms = atoms.repeat()
    if data is None:
        write(filename, atoms, format=format)
    else:
        write(filename, atoms, format=format, data=data)
    if block:
        subprocess.call(command.split() + [filename])
        os.remove(filename)
    else:
        subprocess.Popen(command.split() + [filename])
        subprocess.Popen(['sleep 60; rm {0}'.format(filename)], shell=True)
