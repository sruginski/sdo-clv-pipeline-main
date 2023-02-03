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

df_vels_avg = pd.read_csv(datadir + "average_velocities.csv")
df_vels_avg.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_vels_avg.drop_duplicates()
df_vels_avg.reset_index(drop=True)

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

df_vels_full_avg = df_vels_unw[(np.isnan(df_vels_avg.lo_mu)) & (df_vels_avg.region == 0.0)]
df_vels_full_avg.reset_index(drop=True)
df_vels_full_avg.to_csv(outdir + "full_disk_vels_average.csv", index=False)

df_pixels_full = df_pixels[np.isnan(df_pixels.lo_mu)]
df_pixels_full.reset_index(drop=True)

df_light_full = df_light[np.isnan(df_light.lo_mu)]
df_light_full.reset_index(drop=True)

# get velocities for regions in mu
df_regs = df_vels[(df_vels.region > 0.0) & (~np.isnan(df_vels.lo_mu))]
df_regs.reset_index(drop=True)

df_regs_unw = df_vels_unw[(df_vels_unw.region > 0.0) & (~np.isnan(df_vels_unw.lo_mu))]
df_regs_unw.reset_index(drop=True)

df_regs_avg = df_vels_avg[(df_vels_avg.region > 0.0) & (~np.isnan(df_vels_avg.lo_mu))]
df_regs_avg.reset_index(drop=True)

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# make dfs by mu
plage = df_regs[df_regs.region == 6.0]
network = df_regs[df_regs.region == 5.0]
quiet_sun = df_regs[df_regs.region == 4.0]
red_penumbrae = df_regs[df_regs.region == 3.0]
all_penumbrae = df_regs[df_regs.region == 2.5]
blu_penumbrae = df_regs[df_regs.region == 2.0]
umbrae = df_regs[df_regs.region == 1.0]

# make dfs by mu
plage_unw = df_regs_unw[df_regs_unw.region == 6.0]
network_unw = df_regs_unw[df_regs_unw.region == 5.0]
quiet_sun_unw = df_regs_unw[df_regs_unw.region == 4.0]
red_penumbrae_unw = df_regs_unw[df_regs_unw.region == 3.0]
all_penumbrae_unw = df_regs_unw[df_regs_unw.region == 2.5]
blu_penumbrae_unw = df_regs_unw[df_regs_unw.region == 2.0]
umbrae_unw = df_regs_unw[df_regs_unw.region == 1.0]

# make dfs by mu
plage_avg = df_regs_avg[df_regs_avg.region == 6.0]
network_avg = df_regs_avg[df_regs_avg.region == 5.0]
quiet_sun_avg = df_regs_avg[df_regs_avg.region == 4.0]
red_penumbrae_avg = df_regs_avg[df_regs_avg.region == 3.0]
all_penumbrae_avg = df_regs_avg[df_regs_avg.region == 2.5]
blu_penumbrae_avg = df_regs_avg[df_regs_avg.region == 2.0]
umbrae_avg = df_regs_avg[df_regs_avg.region == 1.0]

# mask rows where all vels are 0.0 (i.e., region isn't present in that annulus)
plage = mask_all_zero_rows(plage)
plage.reset_index(drop=True)
plage.to_csv(outdir + "plage_vels.csv", index=False)
plage_unw = mask_all_zero_rows(plage_unw)
plage_unw.reset_index(drop=True)
plage_unw.to_csv(outdir + "plage_vels_unw.csv", index=False)
plage_avg = mask_all_zero_rows(plage_avg)
plage_avg.reset_index(drop=True)
plage_avg.to_csv(outdir + "plage_vels_avg.csv", index=False)

network = mask_all_zero_rows(network)
network.reset_index(drop=True)
network.to_csv(outdir + "network_vels.csv", index=False)
network_unw = mask_all_zero_rows(network_unw)
network_unw.reset_index(drop=True)
network_unw.to_csv(outdir + "network_vels_unw.csv", index=False)
network_avg = mask_all_zero_rows(network_avg)
network_avg.reset_index(drop=True)
network_avg.to_csv(outdir + "network_vels_avg.csv", index=False)

quiet_sun = mask_all_zero_rows(quiet_sun)
quiet_sun.reset_index(drop=True)
quiet_sun.to_csv(outdir + "quiet_sun_vels.csv", index=False)
quiet_sun_unw = mask_all_zero_rows(quiet_sun_unw)
quiet_sun_unw.reset_index(drop=True)
quiet_sun_unw.to_csv(outdir + "quiet_sun_vels_unw.csv", index=False)
quiet_sun_avg = mask_all_zero_rows(quiet_sun_avg)
quiet_sun_avg.reset_index(drop=True)
quiet_sun_avg.to_csv(outdir + "quiet_sun_vels_avg.csv", index=False)

red_penumbrae = mask_all_zero_rows(red_penumbrae)
red_penumbrae.reset_index(drop=True)
red_penumbrae.to_csv(outdir + "red_penumbrae_vels.csv", index=False)
red_penumbrae_unw = mask_all_zero_rows(red_penumbrae_unw)
red_penumbrae_unw.reset_index(drop=True)
red_penumbrae_unw.to_csv(outdir + "red_penumbrae_vels_unw.csv", index=False)
red_penumbrae_avg = mask_all_zero_rows(red_penumbrae_avg)
red_penumbrae_avg.reset_index(drop=True)
red_penumbrae_avg.to_csv(outdir + "red_penumbrae_vels_avg.csv", index=False)

all_penumbrae = mask_all_zero_rows(all_penumbrae)
all_penumbrae.reset_index(drop=True)
all_penumbrae.to_csv(outdir + "all_penumbrae_vels.csv", index=False)
all_penumbrae_unw = mask_all_zero_rows(all_penumbrae_unw)
all_penumbrae_unw.reset_index(drop=True)
all_penumbrae_unw.to_csv(outdir + "all_penumbrae_vels_unw.csv", index=False)
all_penumbrae_avg = mask_all_zero_rows(all_penumbrae_avg)
all_penumbrae_avg.reset_index(drop=True)
all_penumbrae_avg.to_csv(outdir + "all_penumbrae_vels_avg.csv", index=False)

blu_penumbrae = mask_all_zero_rows(blu_penumbrae)
blu_penumbrae.reset_index(drop=True)
blu_penumbrae.to_csv(outdir + "blu_penumbrae_vels.csv", index=False)
blu_penumbrae_unw = mask_all_zero_rows(blu_penumbrae_unw)
blu_penumbrae_unw.reset_index(drop=True)
blu_penumbrae_unw.to_csv(outdir + "blu_penumbrae_vels_unw.csv", index=False)
blu_penumbrae_avg = mask_all_zero_rows(blu_penumbrae_avg)
blu_penumbrae_avg.reset_index(drop=True)
blu_penumbrae_avg.to_csv(outdir + "blu_penumbrae_vels_avg.csv", index=False)

umbrae = mask_all_zero_rows(umbrae)
umbrae.reset_index(drop=True)
umbrae.to_csv(outdir + "umbrae_vels.csv", index=False)
umbrae_unw = mask_all_zero_rows(umbrae_unw)
umbrae_unw.reset_index(drop=True)
umbrae_unw.to_csv(outdir + "umbrae_vels_unw.csv", index=False)
umbrae_avg = mask_all_zero_rows(umbrae_avg)
umbrae_avg.reset_index(drop=True)
umbrae_avg.to_csv(outdir + "umbrae_vels_avg.csv", index=False)
