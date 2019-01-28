"""LAMMPS has the options to use several internal units (which can be different
from hte ones used in ase).  Mapping is therefore necessary.

See: https://lammps.sandia.gov/doc/units.html
 """

from scipy.constants import gram
from scipy.constants import N_A
from scipy.constants import angstrom
from scipy.constants import eV
from scipy.constants import elementary_charge
from scipy.constants import pico
from scipy.constants import micro
from scipy.constants import atto
from scipy.constants import nano
from scipy.constants import bar
from scipy.constants import centi
from scipy.constants import femto
from scipy.constants import kilo
from scipy.constants import calorie
from scipy.constants import atmosphere
from scipy.constants import erg
from scipy.constants import dyne

from ase import units

PASCAL = 1.0
SECOND = 1.0
METER = 1.0
KELVIN = 1.0
JOULE = 1.0
NEWTON = 1.0
COULOMB = units.C
VOLT = 1.0

ATOMIC_MASS_UNIT = 1.0 / units.kg
MOLE = N_A
POISE = 0.1 * PASCAL * SECOND
STATCOULOMB = 1.0 / 2997924580.0 * COULOMB
STATVOLT = 299.792458 * VOLT

DIM = 3.0

UNITSETS = {}

# !TODO add missing units to 'ASE'
# !TODO add reduce Lennards-Jones units
# !TODO cross check which CODATA version lammps uses
UNITSETS["ASE"] = dict(
    mass=gram / MOLE,
    distance=angstrom,
    time=1.0 / units.second,
    energy=eV,
    velocity=angstrom / (1.0 / units.second),
    force=eV / angstrom,
    pressure=1.0 / units.Pascal,
    charge=elementary_charge,
)

UNITSETS["real"] = dict(
    mass=gram / MOLE,
    distance=angstrom,
    time=femto,
    energy=kilo * calorie / MOLE,
    velocity=angstrom / femto,
    force=kilo * calorie / (MOLE * angstrom),
    torque=kilo * calorie / MOLE,
    temperature=KELVIN,
    pressure=atmosphere,
    dynamic_viscosity=POISE,
    charge=elementary_charge,
    electric_field=1.0 / angstrom,
    density=gram / centi ** DIM,
)

UNITSETS["metal"] = dict(
    mass=gram / MOLE,
    distance=angstrom,
    time=pico,
    energy=eV,
    velocity=angstrom / pico,
    force=eV / angstrom,
    torque=eV,
    temperature=KELVIN,
    pressure=bar,
    dynamic_viscosity=POISE,
    charge=elementary_charge,
    dipole=elementary_charge * angstrom,
    electric_field=1.0 / angstrom,
    density=gram / centi ** DIM,
)

UNITSETS["si"] = dict(
    mass=kilo * gram,
    distance=METER,
    time=SECOND,
    energy=JOULE,
    velocity=METER / SECOND,
    force=NEWTON,
    torque=NEWTON * METER,
    temperature=KELVIN,
    pressure=PASCAL,
    dynamic_viscosity=PASCAL * SECOND,
    charge=COULOMB,
    dipole=COULOMB * METER,
    electric_field=VOLT / METER,
    density=kilo * gram / METER ** DIM,
)

UNITSETS["cgs"] = dict(
    mass=gram,
    distance=centi * METER,
    time=SECOND,
    energy=erg,
    velocity=centi * METER / SECOND,
    force=dyne,
    torque=dyne * centi * METER,
    temperature=KELVIN,
    pressure=dyne / (centi * METER) ** 2,  # or barye = 1.0e-6 bars,
    dynamic_viscosity=POISE,
    charge=STATCOULOMB,  # or esu (4.8032044e-10 is a proton)
    dipole=STATCOULOMB * centi * METER,  # = 10^18 debye,
    electric_field=STATVOLT / (centi * METER),  # or dyne / esu,,
    density=gram / (centi * METER ** DIM),
)

UNITSETS["electron"] = dict(
    mass=ATOMIC_MASS_UNIT,
    distance=units.Bohr,
    time=femto * SECOND,
    energy=units.Hartree / units.J,
    # velocity = Bohr / atomic time units, [1.03275e-15 SECONDs]
    velocity=units.Bohr / (units.AUT / units.second),
    force=units.Hartree / units.Bohr,
    temperature=KELVIN,
    pressure=PASCAL,
    charge=elementary_charge,  # multiple of electron charge (1.0 is a proton)
    dipole=units.Debye,
    electric_field=VOLT / (centi * METER),
)

UNITSETS["micro"] = dict(
    mass=pico * gram,
    distance=micro * METER,
    time=micro * SECOND,
    energy=pico * gram * (micro * METER) ** 2 / (micro * SECOND) ** 2,
    velocity=micro * METER / micro * SECOND,
    force=pico * gram * micro * METER / (micro * SECOND) ** 2,
    torque=pico * gram * (micro * METER) ** 2 / (micro * SECOND) ** 2,
    temperature=KELVIN,
    pressure=pico * gram / (micro * METER * (micro * SECOND) ** 2),
    dynamic_viscosity=pico * gram / (micro * METER * micro * SECOND),
    charge=pico * COULOMB,  # (1.6021765e-7 is a proton),
    dipole=pico * COULOMB * micro * METER,
    electric_field=VOLT / micro * METER,
    density=pico * gram / (micro * METER) ** DIM,
)

UNITSETS["nano"] = dict(
    mass=atto * gram,
    distance=nano * METER,
    time=nano * SECOND,
    energy=atto * gram * (nano * METER) ** 2 / (nano * SECOND) ** 2,
    velocity=nano * METER / (nano * SECOND),
    force=atto * gram * nano * METER / (nano * SECOND) ** 2,
    torque=atto * gram * (nano * METER) ** 2 / (nano * SECOND) ** 2,
    temperature=KELVIN,
    pressure=atto * gram / (nano * METER * (nano * SECOND) ** 2),
    dynamic_viscosity=atto * gram / (nano * METER * nano * SECOND),
    charge=elementary_charge,  # multiple of electron charge (1.0 is a proton)
    dipole=elementary_charge * nano * METER,
    electric_field=VOLT / (nano * METER),
    density=atto * gram / (nano * METER) ** DIM,
)


def convert(value, quantity, fromunits, tounits):
    """Convert units between LAMMPS and ASE.

    :param value: converted value
    :param quantity: mass, distance, time, energy, veloctiy, force, torque,
    temperature, pressure, dynamic_viscosity, charge, dipole,
    electric_field or density
    :param fromunits: ase, metal, real or other (see lammps docs).
    :param tounits: ase, metal, real or other
    :returns: converted value
    :rtype:
    """
    return UNITSETS[fromunits][quantity] / UNITSETS[tounits][quantity] * value
