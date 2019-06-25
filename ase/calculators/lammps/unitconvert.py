"""LAMMPS has the options to use several internal units (which can be different
from the ones used in ase).  Mapping is therefore necessary.

See: https://lammps.sandia.gov/doc/units.html
 """

from .unitconvert_constants import (
    gram_si,
    angstrom_si,
    ev_si,
    e_si,
    picogram_si,
    picocoulomb_si,
    picosecond_si,
    micrometer_si,
    microsecond_si,
    attogram_si,
    nanosecond_si,
    nanometer_si,
    bar_si,
    centimeter_si,
    femtosecond_si,
    kilogram_si,
    gram_per_mole_si,
    ev_per_angstrom_si,
    kcal_per_mole_si,
    angstrom_per_femtosecond_si,
    kcal_per_mole_angstrom_si,
    kelvin_si,
    atmosphere_si,
    volt_per_angstrom_si,
    angstrom_per_picosecond_si,
    poise_si,
    electron_angstrom_si,
    meter_si,
    meter_per_second_si,
    newton_si,
    joule_si,
    pascal_si,
    coulomb_si,
    coulomb_meter_si,
    volt_per_meter_si,
    second_si,
    erg_si,
    centimeter_per_second_si,
    dyne_si,
    dyne_centimeter_si,
    dyne_per_centimetersq_si,
    statcoulomb_si,
    statcoulomb_centimeter_si,
    statvolt_per_centimeter_si,
    amu_si,
    bohr_si,
    hartree_si,
    bohr_per_atu_si,
    hartree_per_bohr_si,
    debye_si,
    volt_per_centimeter_si,
    micrometer_per_microsecond_si,
    picogram_micrometer_per_microsecondsq_si,
    picogram_micrometersq_per_microsecondsq_si,
    picogram_per_micrometer_microsecondsq_si,
    picogram_per_micrometer_microsecond_si,
    picocoulomb_micrometer_si,
    volt_per_micrometer_si,
    nanometer_per_nanosecond_si,
    attogram_nanometer_per_nanosecondsq_si,
    attogram_nanometersq_per_nanosecondsq_si,
    attogram_per_nanometer_nanosecondsq_si,
    attogram_per_nanometer_nanosecond_si,
    electron_nanometer_si,
    volt_si,
)

from ase import units

# !TODO add reduced Lennard-Jones units?

# NOTE: We assume a three-dimensional simulation here!
DIM = 3.0

UNITSETS = {}

UNITSETS["ASE"] = dict(
    mass=1.0 / units.kg,
    distance=1.0 / units.m,
    time=1.0 / units.second,
    energy=1.0 / units.J,
    velocity=units.second / units.m,
    force=units.m / units.J,
    pressure=1.0 / units.Pascal,
    charge=1.0 / units.C,
)

UNITSETS["real"] = dict(
    mass=gram_per_mole_si,
    distance=angstrom_si,
    time=femtosecond_si,
    energy=kcal_per_mole_si,
    velocity=angstrom_per_femtosecond_si,
    force=kcal_per_mole_angstrom_si,
    torque=kcal_per_mole_si,
    temperature=kelvin_si,
    pressure=atmosphere_si,
    dynamic_viscosity=poise_si,
    charge=e_si,
    dipole=electron_angstrom_si,
    electric_field=volt_per_angstrom_si,
    density=gram_si / centimeter_si ** DIM,
)

UNITSETS["metal"] = dict(
    mass=gram_per_mole_si,
    distance=angstrom_si,
    time=picosecond_si,
    energy=ev_si,
    velocity=angstrom_per_picosecond_si,
    force=ev_per_angstrom_si,
    torque=ev_si,
    temperature=kelvin_si,
    pressure=bar_si,
    dynamic_viscosity=poise_si,
    charge=e_si,
    dipole=electron_angstrom_si,
    electric_field=volt_per_angstrom_si,
    density=gram_si / centimeter_si ** DIM,
)

UNITSETS["si"] = dict(
    mass=kilogram_si,
    distance=meter_si,
    time=second_si,
    energy=joule_si,
    velocity=meter_per_second_si,
    force=newton_si,
    torque=joule_si,
    temperature=kelvin_si,
    pressure=pascal_si,
    dynamic_viscosity=pascal_si * second_si,
    charge=coulomb_si,
    dipole=coulomb_meter_si,
    electric_field=volt_per_meter_si,
    density=kilogram_si / meter_si ** DIM,
)

UNITSETS["cgs"] = dict(
    mass=gram_si,
    distance=centimeter_si,
    time=second_si,
    energy=erg_si,
    velocity=centimeter_per_second_si,
    force=dyne_si,
    torque=dyne_centimeter_si,
    temperature=kelvin_si,
    pressure=dyne_per_centimetersq_si,  # or barye = 1.0e-6 bars
    dynamic_viscosity=poise_si,
    charge=statcoulomb_si,  # or esu (4.8032044e-10 is a proton)
    dipole=statcoulomb_centimeter_si,  # = 10^18 debye,
    electric_field=statvolt_per_centimeter_si,  # or dyne / esu
    density=gram_si / (centimeter_si ** DIM),
)

UNITSETS["electron"] = dict(
    mass=amu_si,
    distance=bohr_si,
    time=femtosecond_si,
    energy=hartree_si,
    velocity=bohr_per_atu_si,
    force=hartree_per_bohr_si,
    temperature=kelvin_si,
    pressure=pascal_si,
    charge=e_si,  # multiple of electron charge (1.0 is a proton)
    dipole=debye_si,
    electric_field=volt_per_centimeter_si,
)

UNITSETS["micro"] = dict(
    mass=picogram_si,
    distance=micrometer_si,
    time=microsecond_si,
    energy=picogram_micrometersq_per_microsecondsq_si,
    velocity=micrometer_per_microsecond_si,
    force=picogram_micrometer_per_microsecondsq_si,
    torque=picogram_micrometersq_per_microsecondsq_si,
    temperature=kelvin_si,
    pressure=picogram_per_micrometer_microsecondsq_si,
    dynamic_viscosity=picogram_per_micrometer_microsecond_si,
    charge=picocoulomb_si,  # (1.6021765e-7 is a proton),
    dipole=picocoulomb_micrometer_si,
    electric_field=volt_per_micrometer_si,
    density=picogram_si / (micrometer_si) ** DIM,
)

UNITSETS["nano"] = dict(
    mass=attogram_si,
    distance=nanometer_si,
    time=nanosecond_si,
    energy=attogram_nanometersq_per_nanosecondsq_si,
    velocity=nanometer_per_nanosecond_si,
    force=attogram_nanometer_per_nanosecondsq_si,
    torque=attogram_nanometersq_per_nanosecondsq_si,
    temperature=kelvin_si,
    pressure=attogram_per_nanometer_nanosecondsq_si,
    dynamic_viscosity=attogram_per_nanometer_nanosecond_si,
    charge=e_si,  # multiple of electron charge (1.0 is a proton)
    dipole=electron_nanometer_si,
    electric_field=volt_si / nanometer_si,
    density=attogram_si / nanometer_si ** DIM,
)


def convert(value, quantity, fromunits, tounits):
    """Convert units between LAMMPS and ASE.

    :param value: converted value
    :param quantity: mass, distance, time, energy, velocity, force, torque,
    temperature, pressure, dynamic_viscosity, charge, dipole,
    electric_field or density
    :param fromunits: ASE, metal, real or other (see lammps docs).
    :param tounits: ASE, metal, real or other
    :returns: converted value
    :rtype:
    """
    return UNITSETS[fromunits][quantity] / UNITSETS[tounits][quantity] * value
