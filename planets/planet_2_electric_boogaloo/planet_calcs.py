#%%
"""
Updated version to deal with uncertainty too
"""
#%% Import Libraries and Modules

import astropy.units as u
import astropy.constants as c
import numpy as np

#%% Define Functions

##############################################################################
def spectral_energy_density(wavelength, temperature):
    return (
        2
        * c.h
        * c.c**2
        / wavelength**5
        / (np.exp(c.h * c.c / (wavelength * c.k_B * temperature)) - 1)
    )


##############################################################################
def get_database_value_with_errors(row, key, units=None):

    value = row[key]
    value_min = row[key] - row[key + "_error_min"]
    value_max = row[key] + row[key + "_error_max"]

    if np.isnan(value_min):
        value_min = value
    if np.isnan(value_max):
        value_max = value

    if units is not None:
        value = value * units
        value_min = value_min * units
        value_max = value_max * units

    return [value_min, value, value_max]


##############################################################################
def compute_angular_seperation(database_row):
    
    star_distance_min, star_distance, star_distance_max = (
        get_database_value_with_errors(database_row, "star_distance", u.pc)
    )

    planet_semi_major_axis_min, planet_semi_major_axis, \
        planet_semi_major_axis_max = (
        get_database_value_with_errors(database_row, "semi_major_axis", u.au)
    )

    angular_seperation = np.arctan(
        planet_semi_major_axis / star_distance).to(u.mas)
    
    angular_seperation_min = np.arctan(
        planet_semi_major_axis_min / star_distance_max
    ).to(u.mas)
    
    angular_seperation_max = np.arctan(
        planet_semi_major_axis_max / star_distance_min
    ).to(u.mas)

    return [angular_seperation_min, angular_seperation, angular_seperation_max]


##############################################################################
def mass_to_radius(mass, density=1.64 * u.g / (u.cm**3)):

    volume = mass / density
    radius = (3 * volume / (4 * np.pi)) ** (1 / 3)

    return radius


##############################################################################
def compute_contrast(database_row, wavel, albedo=0.0,
                      density=1.64 * u.g / u.cm**3):
    
    # # set up function for spectral energy density
    # spectral_energy_density = spectral_energy_density
    
    star_temp = database_row["star_teff"] * u.K  # no errors
    star_radius = database_row["star_radius"] * u.Rsun  # no errors

    if np.isnan(database_row["mass"]):
        # use mass_sin_i, the minimum mass of the planet due to 
        # inclination effect
        planet_mass_min, planet_mass, planet_mass_max = (
            get_database_value_with_errors(
                database_row, "mass_sin_i", units=u.Mjup
            )
        )
        uses_sin_i = True

    else:
        planet_mass_min, planet_mass, planet_mass_max = (
            get_database_value_with_errors(
                database_row, "mass", units=u.Mjup
            )
        )
        uses_sin_i = False

    # planet_mass_min, planet_mass, planet_mass_max = planet_calcs.get_database_value_with_errors(database_row, "mass_sin_i", units=u.Mjup)
    if np.isnan(planet_mass):
        return [np.nan * u.dimensionless_unscaled] * 3, False

    # infer radius
    planet_radius_min = mass_to_radius(planet_mass_min, density)
    planet_radius = mass_to_radius(planet_mass, density)
    planet_radius_max = mass_to_radius(planet_mass_max, density)

    star_bolometric_luminosity = 4 * np.pi * star_radius**2 \
        * c.sigma_sb * star_temp**4

    # now eq temp
    semi_major_axis_min, semi_major_axis, semi_major_axis_max = (
        get_database_value_with_errors(
            database_row, "semi_major_axis", units=u.AU
        )
    )

    planet_T_eq = (
        star_bolometric_luminosity
        * (1 - albedo)
        / (16 * np.pi * c.sigma_sb * semi_major_axis**2)
    ) ** (1 / 4)
    planet_T_eq_min = (
        star_bolometric_luminosity
        * (1 - albedo)
        / (16 * np.pi * c.sigma_sb * semi_major_axis_max**2)
    ) ** (1 / 4)
    planet_T_eq_max = (
        star_bolometric_luminosity
        * (1 - albedo)
        / (16 * np.pi * c.sigma_sb * semi_major_axis_min**2)
    ) ** (1 / 4)

    star_energy_density = spectral_energy_density(wavel, star_temp)
    planet_energy_density = spectral_energy_density(wavel, planet_T_eq)
    planet_energy_density_min = spectral_energy_density(wavel, planet_T_eq_min)
    planet_energy_density_max = spectral_energy_density(wavel, planet_T_eq_max)

    ratio = (planet_energy_density * np.pi * planet_radius**2) / (
        star_energy_density * np.pi * star_radius**2
    )
    ratio_min = (planet_energy_density_min * np.pi * planet_radius_min**2) / (
        star_energy_density * np.pi * star_radius**2
    )
    ratio_max = (planet_energy_density_max * np.pi * planet_radius_max**2) / (
        star_energy_density * np.pi * star_radius**2
    )

    return [ratio_min.to(""), ratio.to(""), ratio_max.to("")], uses_sin_i
