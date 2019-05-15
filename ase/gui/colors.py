# -*- coding: utf-8 -*-
"""colors.py - select how to color the atoms in the GUI."""
from __future__ import unicode_literals
from ase.gui.i18n import _

import numpy as np

import ase.gui.ui as ui
from ase.gui.utils import get_magmoms


class ColorWindow:
    """A window for selecting how to color the atoms."""
    def __init__(self, gui):
        self.win = ui.Window(_('Colors'))
        self.gui = gui
        self.win.add(ui.Label(_('Choose how the atoms are colored:')))
        values = ['jmol', 'tag', 'force', 'velocity',
                  'initial charge', 'magmom', 'neighbors']
        labels = [_('By atomic number, default "jmol" colors'),
                  _('By tag'),
                  _('By force'),
                  _('By velocity'),
                  _('By initial charge'),
                  _('By magnetic moment'),
                  _('By number of neighbors'),]

        self.radio = ui.RadioButtons(labels, values, self.toggle,
                                     vertical=True)
        self.radio.value = gui.colormode
        self.win.add(self.radio)
        self.activate()
        self.label = ui.Label()
        self.win.add(self.label)
        if hasattr(self, 'mnmx'):
            self.win.add(self.mnmx)

    def change_mnmx(self):
        ##print('change min/max value')
        mn, mx = self.mnmx[1].value, self.mnmx[3].value
        colorscale, _, _ = self.gui.colormode_data
        self.gui.colormode_data = colorscale, mn, mx
        self.gui.draw()

    def activate(self):
        images = self.gui.images
        atoms = self.gui.atoms
        radio = self.radio
        radio['tag'].active = atoms.has('tags')

        # XXX not sure how to deal with some images having forces,
        # and other images not.  Same goes for below quantities
        F = images.get_forces(atoms)
        radio['force'].active = F is not None
        radio['velocity'].active = atoms.has('momenta')
        radio['initial charge'].active = atoms.has('initial_charges')
        radio['magmom'].active = get_magmoms(atoms).any()
        radio['neighbors'].active = True

    def toggle(self, value):
        self.gui.colormode = value
        if value == 'jmol' or value == 'neighbors':
            if hasattr(self, 'mnmx'):
                "delete the min max fields by creating a new window"
                del self.mnmx
                self.win.close()
                self.__init__(self.gui)
            text = ''
        else:
            scalars = np.array([self.gui.get_color_scalars(i)
                                for i in range(len(self.gui.images))])
            mn = np.nanmin(scalars)
            mx = np.nanmax(scalars)
            colorscale = ['#{0:02X}80{0:02X}'.format(red)
                          for red in range(0, 260, 10)]
            self.gui.colormode_data = colorscale, mn, mx

            unit = {'tag': '',
                    'force': 'eV/Ang',
                    'velocity': '??',
                    'charge': '|e|',
                    'initial charge': '|e|',
                    u'magmom': 'μB'}[value]
            text = '[{0},{1}]: [{2:.6f},{3:.6f}] {4}'.format(
                _('Green'), _('Yellow'), mn, mx, unit)
            text = ''

            rng = mx - mn  # XXX what are optimal allowed range and steps ?
            self.mnmx = [_('min:'),
                         ui.SpinBox(mn, mn - 10 * rng, mx + rng, rng / 10.,
                                    self.change_mnmx),
                         _('max:'),
                         ui.SpinBox(mx, mn - 10 * rng, mx + rng, rng / 10.,
                                    self.change_mnmx),
                         _(unit)
            ]
            self.win.close()
            self.__init__(self.gui)

        self.label.text = text
        self.radio.value = value
        self.gui.draw()
        return text  # for testing

    def notify_atoms_changed(self):
        "Called by gui object when the atoms have changed."
        self.activate()
        mode = self.gui.colormode
        if not self.radio[mode].active:
            mode = 'jmol'
        self.toggle(mode)
