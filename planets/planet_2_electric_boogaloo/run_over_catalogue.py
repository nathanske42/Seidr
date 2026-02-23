#%%
import numpy as np
import time
import pyvo

import planet_calcs

import astropy.units as u
import astropy.constants as c

import matplotlib.pyplot as plt


#%% Import Planet Catalogue

# docs: http://voparis-tap-planeto.obspm.fr/__system__/dc_tables/show/tableinfo/exoplanet.epn_core
service = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")

query = "SELECT * FROM exoplanet.epn_core"

start = time.time()
results = service.search(query)
end = time.time()

print(f"Query took {end - start} seconds")

print(len(results))
print(len(results[0]))

np.save("myresults", results.to_table())

pandas_table = results.to_table().to_pandas()

#%%

# lambda
wavel = 1.65 * u.micron

## Initialize lists to store results
seperations = []
sep_errors = []
contrasts = []
contrast_errors = []
sin_i_mask = []


for index, row in pandas_table.iterrows():

    sep = planet_calcs.compute_angular_seperation(row)
    contrast, uses_sin_i = planet_calcs.compute_contrast(row, wavel)

    sep_err = [sep[1].value - sep[0].value, sep[2].value - sep[1].value]
    if (np.array(sep_err) < 0).sum() > 0:
        seperations.append(np.nan)
        sep_errors.append([np.nan, np.nan])
    else:
        seperations.append(sep[1].value)
        sep_errors.append(sep_err)

    c_err = [
        contrast[1].value - contrast[0].value,
        contrast[2].value - contrast[1].value,
    ]

    if (np.array(c_err) < 0).sum() > 0:
        contrasts.append(np.nan)
        contrast_errors.append([np.nan, np.nan])
    else:
        contrasts.append(contrast[1].value)
        contrast_errors.append(c_err)

    sin_i_mask.append(uses_sin_i)


## Add to pandas table
pandas_table["contrast"] = contrasts
pandas_table["seperation"] = seperations


#%% Convert to Numpy Arrays

sep_errors = np.array(sep_errors).T
contrast_errors = np.array(contrast_errors).T
seperations = np.array(seperations)
contrasts = np.array(contrasts)
sin_i_mask = np.array(sin_i_mask)

sep_errors.shape, sin_i_mask.shape


#%%

np.isnan(sep_errors).sum(axis=0), np.isnan(contrast_errors).sum(axis=0)
has_errors = np.logical_and(
    np.isnan(sep_errors).sum(axis=0) == 0, 
    np.isnan(contrast_errors).sum(axis=0) == 0
)

#%%
sep_errors[:, has_errors]
print("Number of planets with valid separation and contrast errors:", has_errors.sum())

#%% Check Declination Angles

pandas_table["dec"]

is_good_dec = np.array([pandas_table["dec"] < 20]).flatten()
is_good_dec

#%% Plot Contrast vs Separation

# plt.scatter(seperations, contrasts, s=5)
# plt.errorbar(seperations[sin_i_mask], contrasts[sin_i_mask], xerr=sep_errors[sin_i_mask], yerr=contrast_errors[sin_i_mask], fmt="o", markersize=5,)

# sin_i_mask = np.logical_not(sin_i_mask)
# plt.errorbar(seperations[sin_i_mask], contrasts[sin_i_mask], xerr=sep_errors[sin_i_mask], yerr=contrast_errors[sin_i_mask], fmt="o", markersize=5,)
plt.errorbar(
    seperations[is_good_dec],
    contrasts[is_good_dec],
    xerr=sep_errors[:, is_good_dec],
    yerr=contrast_errors[:, is_good_dec],
    fmt="o",
    markersize=3,
    label="dec < 20°",
)

is_good_dec = np.logical_not(is_good_dec)

plt.errorbar(
    seperations[is_good_dec],
    contrasts[is_good_dec],
    xerr=sep_errors[:, is_good_dec],
    yerr=contrast_errors[:, is_good_dec],
    fmt="o",
    markersize=3,
    label="dec > 20°",
    alpha=0.2,
)

is_good_dec = np.logical_not(is_good_dec)


longest_baseline = 130 * u.m
contrast_ylimit = 10**-8

plt.xlim(0.9, 12)
plt.yscale("log")
# plt.xscale("log")
plt.ylim(contrast_ylimit, 10**-3)
plt.axvline(
    ((wavel / (2 * longest_baseline)) * u.rad).to(u.mas).value,
    color="black",
    linestyle="--",
    label="0.5 λ/longest baseline",
)
plt.legend()
plt.xlabel("Angular separation (milliarcseconds)")
plt.ylabel("Contrast of thermal emission in H band")
plt.savefig("contrast_vs_seperation.pdf")


# %%
pandas_table["log10contrast"] = np.log10(pandas_table["contrast"])


#%% Extract Achievable Planets

achievable = pandas_table[pandas_table["contrast"] > 5 * 10**-6]
achievable = achievable[
    achievable["seperation"]
    > ((wavel / (2 * longest_baseline)) * u.rad).to(u.mas).value
]
achievable = achievable[achievable["dec"] < 20]

table = achievable[["target_name", "log10contrast", "seperation", "dec"]]
print(
    table.sort_values("log10contrast", ascending=False).to_latex(
        index=False, float_format="{:.2f}".format
    )
)

# %%
