"""Build ASE's web-page."""

import os
import shutil
import subprocess
import sys
from pathlib import Path


cmds = """\
python3 -m venv venv
. venv/bin/activate
pip install sphinx-rtd-theme pillow
git clone http://gitlab.com/ase_koopmans/ase_koopmans.git
cd ase_koopmans
git checkout {branch}
pip install .
python setup.py sdist
cd doc
make
mv build/html ase_koopmans-web-page"""


def build(branch='master'):
    root = Path(f'/tmp/ase_koopmans-docs-{branch}')
    if root.is_dir():
        sys.exit('Locked')
    root.mkdir()
    os.chdir(root)
    cmds2 = ' && '.join(cmds.format(branch=branch).splitlines())
    p = subprocess.run(cmds2, shell=True)
    if p.returncode == 0:
        status = 'ok'
    else:
        print('FAILED!', file=sys.stdout)
        status = 'error'
    f = root.with_name(f'ase_koopmans-docs-{branch}-{status}')
    if f.is_dir():
        shutil.rmtree(f)
    root.rename(f)
    return status


def build_both():
    assert build('master') == 'ok'
    assert build('web-page') == 'ok'
    tar = next(
        Path('/tmp/ase_koopmans-docs-master-ok/ase_koopmans/dist/').glob('ase_koopmans-*.tar.gz'))
    master = Path('/tmp/ase_koopmans-docs-master-ok/ase_koopmans/doc/ase_koopmans-web-page')
    webpage = Path('/tmp/ase_koopmans-docs-web-page-ok/ase_koopmans/doc/ase_koopmans-web-page')
    home = Path.home() / 'web-pages'
    cmds = ' && '.join(
        [f'cp -rp {master} {webpage}/dev',
         f'cp {tar} {webpage}',
         f'cp {tar} {webpage}/dev',
         f'find {webpage} -name install.html | '
         f'xargs sed -i s/snapshot.tar.gz/{tar.name}/g',
         f'cd {webpage.parent}',
         f'tar -czf ase_koopmans-web-page.tar.gz ase_koopmans-web-page',
         f'cp ase_koopmans-web-page.tar.gz {home}'])
    subprocess.run(cmds, shell=True, check=True)


if __name__ == '__main__':
    # build()
    build_both()
