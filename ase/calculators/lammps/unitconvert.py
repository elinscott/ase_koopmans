"""LAMMPS has the options to use several internal units (which can be different
from hte ones used in ase).  Mapping is therefore necessary
 """

from scipy.constants import gram
from scipy.constants import N_A
from scipy.constants import angstrom
from scipy.constants import eV
from scipy.constants import elementary_charge
from scipy.constants import pico
from scipy.constants import bar
from scipy.constants import centi
from scipy.constants import femto
from scipy.constants import kilo
from scipy.constants import calorie
from scipy.constants import atmosphere

from ase import units
mole = N_A

unitSets = {}

# !TODO add missing units
unitSets['ASE'] = dict(
    mass=gram/mole,
    distance=angstrom,
    time=1./units.second,
    energy=eV,
    velocity=angstrom/(1./units.second),
    force=eV/angstrom,
    pressure=1./units.Pascal,
    stress=1./units.Pascal,
    charge=elementary_charge
)

unitSets['metal'] = dict(
    mass=gram/mole,
    distance=angstrom,
    time=pico,
    energy=eV,
    velocity=angstrom/pico,
    force=eV/angstrom,
    torque=eV,
    temperature=1.,
    pressure=bar,
    stress=bar,
    charge=elementary_charge,
    dipole=elementary_charge*angstrom,
    electricField=1./angstrom,
    density=gram/centi**3.
)

unitSets['real'] = dict(
    mass=gram/mole,
    distance=angstrom,
    time=femto,
    energy=kilo*calorie/mole,
    velocity=angstrom/femto,
    force=kilo*calorie/(mole*angstrom),
    torque=kilo*calorie/mole,
    temperature=1.,
    pressure=atmosphere,
    stress=atmosphere,
    charge=elementary_charge,
    electricField=1./angstrom,
    density=gram/centi**3.
)


def convert(value, quantity, fromunits, tounits):
    """Convert units between LAMMPS and ASE.

    :param value: converted value
    :param quantity: mass, distance, time, energy, veloctiy, force, torque,
    temperature, pressure, stress, charge, dipole,
    electricField or density
    :param fromunits: ase, metal or real
    :param tounits: ase, metal or real
    :returns: converted value
    :rtype:
    """
    return unitSets[fromunits][quantity]/unitSets[tounits][quantity] * value
