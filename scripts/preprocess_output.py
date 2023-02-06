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

# read in and sort by mjd
df_all = pd.read_csv(datadir + "region_output.csv")
df_all.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_all.drop_duplicates()
df_all.reset_index(drop=True, inplace=True)

# get full disk only
df_full_disk = df_all[(np.isnan(df_all.lo_mu)) & np.isnan(df_all.region)]
df_full_disk.reset_index(drop=True, inplace=True)
df_full_disk.to_csv(outdir + "full_disk.csv", index=False)

"""
# find dates with extreme outliers
v_conv_rolling_avg = df_full_disk.v_conv.rolling(30).mean()
v_conv_rolling_std = df_full_disk.v_conv.rolling(30).std()

# find distance from rolling avergge
dist = np.abs(df_full_disk.v_conv - v_conv_rolling_avg)
idx = dist[dist > 2.0 * v_conv_rolling_std].index
"""

# make dfs by mu
plage = df_all[df_all.region == 6.0]
network = df_all[df_all.region == 5.0]
quiet_sun = df_all[df_all.region == 4.0]
red_penumbrae = df_all[df_all.region == 3.0]
all_penumbrae = df_all[df_all.region == 2.5]
blu_penumbrae = df_all[df_all.region == 2.0]
umbrae = df_all[df_all.region == 1.0]

# mask rows where all vels are 0.0 (i.e., region isn't present in that annulus)
plage = mask_all_zero_rows(plage)
plage.reset_index(drop=True, inplace=True)
plage.to_csv(outdir + "plage.csv", index=False)

network = mask_all_zero_rows(network)
network.reset_index(drop=True, inplace=True)
network.to_csv(outdir + "network.csv", index=False)

quiet_sun = mask_all_zero_rows(quiet_sun)
quiet_sun.reset_index(drop=True, inplace=True)
quiet_sun.to_csv(outdir + "quiet_sun.csv", index=False)

red_penumbrae = mask_all_zero_rows(red_penumbrae)
red_penumbrae.reset_index(drop=True, inplace=True)
red_penumbrae.to_csv(outdir + "red_penumbrae.csv", index=False)

all_penumbrae = mask_all_zero_rows(all_penumbrae)
all_penumbrae.reset_index(drop=True, inplace=True)
all_penumbrae.to_csv(outdir + "penumbrae.csv", index=False)

blu_penumbrae = mask_all_zero_rows(blu_penumbrae)
blu_penumbrae.reset_index(drop=True, inplace=True)
blu_penumbrae.to_csv(outdir + "blu_penumbrae.csv", index=False)

umbrae = mask_all_zero_rows(umbrae)
umbrae.reset_index(drop=True, inplace=True)
umbrae.to_csv(outdir + "umbrae.csv", index=False)
