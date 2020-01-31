import argparse
import os
from functools import partial
import unittest

import pytest


if not os.environ.get('DISPLAY'):
    raise unittest.SkipTest('No display')

try:
    import tkinter  # noqa
except ImportError:
    raise unittest.SkipTest('tkinter not available')


import numpy as np

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.build import molecule
import ase.gui.ui as ui
from ase.gui.i18n import _
from ase.gui.gui import GUI
from ase.gui.save import save_dialog


class Error:
    """Fake window for testing puposes."""
    has_been_called = False

    def __call__(self, title, text=None):
        self.text = text or title
        self.has_been_called = True

    def called(self, text=None):
        """Check that an oops-window was opened with correct title."""
        if not self.has_been_called:
            return False

        self.has_been_called = False  # ready for next call

        return text is None or text == self.text


ui.error = Error()


@pytest.fixture
def gui():
    gui = GUI()
    yield gui
    gui.exit()


def test_nanotube(gui):
    nt = gui.nanotube_window()
    nt.apply()
    nt.element[1].value = '?'
    nt.apply()
    assert ui.error.called(
        _('You have not (yet) specified a consistent set of parameters.'))

    nt.element[1].value = 'C'
    nt.ok()
    assert len(gui.images[0]) == 20


def test_nanopartickle(gui):
    n = gui.nanoparticle_window()
    n.element.symbol = 'Cu'
    n.apply()
    n.set_structure_data()
    assert len(gui.images[0]) == 675
    n.method.value = 'wulff'
    n.update_gui_method()
    n.apply()


def test_color(gui):
    a = Atoms('C10', magmoms=np.linspace(1, -1, 10))
    a.positions[:] = np.linspace(0, 9, 10)[:, None]
    a.calc = SinglePointCalculator(a, forces=a.positions)
    che = np.linspace(100, 110, 10)
    mask = [0] * 10
    mask[5] = 1
    a.set_array('corehole_energies', np.ma.array(che, mask=mask))
    gui.new_atoms(a)
    c = gui.colors_window()
    c.toggle('force')
    c.toggle('magmom')
    activebuttons = [button.active for button in c.radio.buttons]
    assert activebuttons == [1, 0, 1, 0, 0, 1, 1, 1], activebuttons
    c.toggle('corehole_energies')
    c.change_mnmx(101, 120)


def test_settings(gui):
    gui.new_atoms(molecule('H2O'))
    s = gui.settings()
    s.scale.value = 1.9
    s.scale_radii()


def test_rotate(gui):
    gui.window['toggle-show-bonds'] = True
    gui.new_atoms(molecule('H2O'))
    gui.rotate_window()


def test_open_and_save(gui):
    mol = molecule('H2O')
    for i in range(3):
        mol.write('h2o.json')
    gui.open(filename='h2o.json')
    save_dialog(gui, 'h2o.cif@-1')


def test_fracocc(gui):
    from ase.test.fio.cif import content
    with open('./fracocc.cif', 'w') as f:
        f.write(content)
    gui.open(filename='fracocc.cif')


def window():

    def hello(event=None):
        print('hello', event)

    menu = [('Hi', [ui.MenuItem('_Hello', hello, 'Ctrl+H')]),
            ('Hell_o', [ui.MenuItem('ABC', hello, choices='ABC')])]
    win = ui.MainWindow('Test', menu=menu)

    win.add(ui.Label('Hello'))
    win.add(ui.Button('Hello', hello))

    r = ui.Rows([ui.Label(x * 7) for x in 'abcd'])
    win.add(r)
    r.add('11111\n2222\n333\n44\n5')

    def abc(x):
        print(x, r.rows)

    cb = ui.ComboBox(['Aa', 'Bb', 'Cc'], callback=abc)
    win.add(cb)

    rb = ui.RadioButtons(['A', 'B', 'C'], 'ABC', abc)
    win.add(rb)

    b = ui.CheckButton('Hello')

    def hi():
        print(b.value, rb.value, cb.value)
        del r[2]
        r.add('-------------')

    win.add([b, ui.Button('Hi', hi)])

    return win


def runcallbacks(win):
    win.things[1].callback()
    win.things[1].callback()
    win.close()


def test_callbacks():
    win = window()
    win.win.after_idle(runcallbacks)
