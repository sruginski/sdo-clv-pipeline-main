import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

from sdo_pypline.paths import root

def mask_all_zero_rows(df, return_idx=False):
    idx = (df.v_hat == 0.0) & (df.v_phot == 0.0) & (df.v_conv == 0.0) & (df.v_quiet == 0.0)
    if return_idx:
        return df[~idx], ~(idx.values)
    else:
        return df[~idx]
    return None

def mask_zero_v_conv(df):
    idx = (df.v_hat == df.v_quiet)
    return df[~idx]

# sort out paths
datadir = str(root / "data") + "/"

# make directory for processed output
if not os.path.isdir(datadir + "processed/"):
    os.mkdir(datadir + "processed/")

outdir = datadir + "processed/"

# read in the data and sort by mjd
df_vels = pd.read_csv(datadir + "velocities.csv")
df_vels.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_vels.drop_duplicates()
df_vels.reset_index(drop=True)

df_vels_unw = pd.read_csv(datadir + "unweighted_velocities.csv")
df_vels_unw.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_vels_unw.drop_duplicates()
df_vels_unw.reset_index(drop=True)

df_pixels = pd.read_csv(datadir + "pixel_stats.csv")
df_pixels.sort_values(by=["mjd", "lo_mu"], inplace=True)
df_pixels.drop_duplicates()
df_pixels.reset_index(drop=True)

df_light = pd.read_csv(datadir + "light_stats.csv")
df_light.sort_values(by=["mjd", "lo_mu"], inplace=True)
df_light.drop_duplicates()
df_light.reset_index(drop=True)

df_intensities = pd.read_csv(datadir + "intensities.csv")
df_intensities.sort_values(by=["mjd"], inplace=True)
df_intensities.drop_duplicates()
df_intensities.reset_index(drop=True)

# get numbers computed for full disk
df_vels_full = df_vels[(np.isnan(df_vels.lo_mu)) & (df_vels.region == 0.0)]
df_vels_full.reset_index(drop=True)
df_vels_full.to_csv(outdir + "full_disk_vels.csv", index=False)

df_vels_full_unw = df_vels_unw[(np.isnan(df_vels_unw.lo_mu)) & (df_vels_unw.region == 0.0)]
df_vels_full_unw.reset_index(drop=True)
df_vels_full_unw.to_csv(outdir + "full_disk_vels_unweighted.csv", index=False)

df_pixels_full = df_pixels[np.isnan(df_pixels.lo_mu)]
df_pixels_full.reset_index(drop=True)

df_light_full = df_light[np.isnan(df_light.lo_mu)]
df_light_full.reset_index(drop=True)

# get velocities for regions in mu
df_regs = df_vels[(df_vels.region > 0.0) & (~np.isnan(df_vels.lo_mu))]
df_regs.reset_index(drop=True)

# get velocities for regions in mu
df_regs_unw = df_vels_unw[(df_vels_unw.region > 0.0) & (~np.isnan(df_vels_unw.lo_mu))]
df_regs_unw.reset_index(drop=True)

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# make dfs by mu
plage = df_regs[df_regs.region == 6.0]
network = df_regs[df_regs.region == 5.0]
quiet_sun = df_regs[df_regs.region == 4.0]
red_penumbrae = df_regs[df_regs.region == 3.0]
blu_penumbrae = df_regs[df_regs.region == 2.0]
umbrae = df_regs[df_regs.region == 1.0]

# make dfs by mu
plage_unw = df_regs_unw[df_regs_unw.region == 6.0]
network_unw = df_regs_unw[df_regs_unw.region == 5.0]
quiet_sun_unw = df_regs_unw[df_regs_unw.region == 4.0]
red_penumbrae_unw = df_regs_unw[df_regs_unw.region == 3.0]
blu_penumbrae_unw = df_regs_unw[df_regs_unw.region == 2.0]
umbrae_unw = df_regs_unw[df_regs_unw.region == 1.0]

# mask rows where all vels are 0.0 (i.e., region isn't present in that annulus)
plage = mask_all_zero_rows(plage)
plage.reset_index(drop=True)
plage.to_csv(outdir + "plage_vels.csv", index=False)
plage_unw = mask_all_zero_rows(plage)
plage_unw.reset_index(drop=True)
plage_unw.to_csv(outdir + "plage_vels_unw.csv", index=False)

network = mask_all_zero_rows(network)
network.reset_index(drop=True)
network.to_csv(outdir + "network_vels.csv", index=False)
network_unw = mask_all_zero_rows(network)
network_unw.reset_index(drop=True)
network_unw.to_csv(outdir + "network_vels_unw.csv", index=False)

quiet_sun = mask_all_zero_rows(quiet_sun)
quiet_sun.reset_index(drop=True)
quiet_sun.to_csv(outdir + "quiet_sun_vels.csv", index=False)
quiet_sun_unw = mask_all_zero_rows(quiet_sun)
quiet_sun_unw.reset_index(drop=True)
quiet_sun_unw.to_csv(outdir + "quiet_sun_vels_unw.csv", index=False)

red_penumbrae = mask_all_zero_rows(red_penumbrae)
red_penumbrae.reset_index(drop=True)
red_penumbrae.to_csv(outdir + "red_penumbrae_vels.csv", index=False)
red_penumbrae_unw = mask_all_zero_rows(red_penumbrae)
red_penumbrae_unw.reset_index(drop=True)
red_penumbrae_unw.to_csv(outdir + "red_penumbrae_vels_unw.csv", index=False)

blu_penumbrae = mask_all_zero_rows(blu_penumbrae)
blu_penumbrae.reset_index(drop=True)
blu_penumbrae.to_csv(outdir + "blu_penumbrae_vels.csv", index=False)
blu_penumbrae_unw = mask_all_zero_rows(blu_penumbrae)
blu_penumbrae_unw.reset_index(drop=True)
blu_penumbrae_unw.to_csv(outdir + "blu_penumbrae_vels_unw.csv", index=False)

umbrae = mask_all_zero_rows(umbrae)
umbrae.reset_index(drop=True)
umbrae.to_csv(outdir + "umbrae_vels.csv", index=False)
umbrae_unw = mask_all_zero_rows(umbrae)
umbrae_unw.reset_index(drop=True)
umbrae_unw.to_csv(outdir + "umbrae_vels_unw.csv", index=False)

# combine red and blue_penumbrae into total penumbrae
pen_light = df_light[~np.isnan(df_light.lo_mu)]
pen_light.reset_index(drop=True)

# make data frame to hold vels
all_penumbrae = pd.DataFrame(columns = red_penumbrae.columns.values)
all_penumbrae_unw = pd.DataFrame(columns = red_penumbrae_unw.columns.values)

# loop over dates or something
for i, mjd in enumerate(np.unique(pen_light.mjd)):
    for j, mu in enumerate(np.unique(pen_light.lo_mu)):
        # get light fractions
        row_light = pen_light[(pen_light.mjd == mjd) & (pen_light.lo_mu == mu)]
        red_light = row_light.red_pen_frac.values[0]
        blu_light = row_light.blu_pen_frac.values[0]
        tot_light = red_light + blu_light

        # get velocities
        red_vels = red_penumbrae[(red_penumbrae.mjd == mjd) & (red_penumbrae.lo_mu == mu)]
        blu_vels = blu_penumbrae[(blu_penumbrae.mjd == mjd) & (blu_penumbrae.lo_mu == mu)]
        red_vels_unw = red_penumbrae_unw[(red_penumbrae_unw.mjd == mjd) & (red_penumbrae_unw.lo_mu == mu)]
        blu_vels_unw = blu_penumbrae_unw[(blu_penumbrae_unw.mjd == mjd) & (blu_penumbrae_unw.lo_mu == mu)]

        # calculate
        if (len(red_vels) == 0) & (len(blu_vels) == 0):
            continue
        elif (len(red_vels) == 0) & (len(blu_vels) != 0):
            v_hat = ((blu_vels.v_hat.values * blu_light) / tot_light)[0]
            v_phot = ((blu_vels.v_phot.values * blu_light) / tot_light)[0]
            v_quiet = ((blu_vels.v_quiet.values * blu_light) / tot_light)[0]
            v_conv = ((blu_vels.v_conv.values * blu_light) / tot_light)[0]

            v_hat_unw = blu_vels_unw.v_hat.values
            v_phot_unw = blu_vels_unw.v_phot.values
            v_quiet_unw = blu_vels_unw.v_quiet.values
            v_conv_unw = blu_vels_unw.v_conv.values
        elif (len(red_vels) != 0) & (len(blu_vels) == 0):
            v_hat = ((red_vels.v_hat.values * red_light) / tot_light)[0]
            v_phot = ((red_vels.v_phot.values * red_light) / tot_light)[0]
            v_quiet = ((red_vels.v_quiet.values * red_light) / tot_light)[0]
            v_conv = ((red_vels.v_conv.values * red_light) / tot_light)[0]

            v_hat_unw = red_vels_unw.v_hat.values
            v_phot_unw = red_vels_unw.v_phot.values
            v_quiet_unw = red_vels_unw.v_quiet.values
            v_conv_unw = red_vels_unw.v_conv.values
        else:
            v_hat = ((red_vels.v_hat.values * red_light + blu_vels.v_hat.values * blu_light) / tot_light)[0]
            v_phot = ((red_vels.v_phot.values * red_light + blu_vels.v_phot.values * blu_light) / tot_light)[0]
            v_quiet = ((red_vels.v_quiet.values * red_light + blu_vels.v_quiet.values * blu_light) / tot_light)[0]
            v_conv = ((red_vels.v_conv.values * red_light + blu_vels.v_conv.values * blu_light) / tot_light)[0]

            v_hat_unw = (red_vels_unw.v_hat.values + blu_vels_unw.v_hat.values)/2.0
            v_phot_unw = (red_vels_unw.v_phot.values + blu_vels_unw.v_phot.values)/2.0
            v_quiet_unw = (red_vels_unw.v_quiet.values + blu_vels_unw.v_quiet.values)/2.0
            v_conv_unw = (red_vels_unw.v_conv.values + blu_vels_unw.v_conv.values)/2.0


        # append to dataframe
        row = pd.Series({"mjd": mjd, "region": 2.5, "lo_mu": mu, "hi_mu": mu + 0.1,
                         "v_hat": v_hat, "v_phot": v_phot, "v_quiet": v_quiet, "v_conv": v_conv})
        row_unw = pd.Series({"mjd": mjd, "region": 2.5, "lo_mu": mu, "hi_mu": mu + 0.1,
                             "v_hat": v_hat_unw, "v_phot": v_phot_unw,
                             "v_quiet": v_quiet_unw, "v_conv": v_conv_unw})
        all_penumbrae = pd.concat([all_penumbrae, pd.DataFrame([row], columns=row.index)]).reset_index(drop=True)
        all_penumbrae_unw = pd.concat([all_penumbrae_unw, pd.DataFrame([row], columns=row.index)]).reset_index(drop=True)

# write to disk
all_penumbrae.reset_index(drop=True)
all_penumbrae.to_csv(outdir + "penumbrae_vels.csv", index=False)

all_penumbrae_unw.reset_index(drop=True)
all_penumbrae_unw.to_csv(outdir + "penumbrae_vels_unw.csv", index=False)

# concatenate to full data frame
df_vels_new = pd.concat([all_penumbrae, df_vels]).reset_index(drop=True)
df_vels_new.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_vels_new.drop_duplicates()
df_vels_new.reset_index(drop=True)
df_vels_new.to_csv(outdir + "velocities.csv", index=False)

# concatenate to full data frame
df_vels_unw_new = pd.concat([all_penumbrae_unw, df_vels_unw]).reset_index(drop=True)
df_vels_unw_new.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_vels_unw_new.drop_duplicates()
df_vels_unw_new.reset_index(drop=True)
df_vels_unw_new.to_csv(outdir + "unweighted_velocities.csv", index=False)
