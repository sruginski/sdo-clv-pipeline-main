import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

from sdo_pypline.paths import root

def mask_all_zero_rows(df):
    idx = (df.mag_avg == 0.0) & (df.mag_std == 0.0) & (df.mag_unsigned == 0.0)
    return df[~idx]

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# read in the data and sort by mjd
df_all = pd.read_csv(datadir + "mag_stats.csv")
df_all.sort_values(by=["mjd"], inplace=True)

# get full disk velocities
df_full = df_all[(np.isnan(df_all.lo_mu)) & (df_all.region == 0.0)]

# get velocites for full mu annuli
df_mu = df_all[(~np.isnan(df_all.lo_mu)) & (df_all.region == 0.0)]

# get velocities for regions in mu
df_regs = df_all[(df_all.region > 0.0) & (~np.isnan(df_all.lo_mu))]

# get velocities for regions on fulldisk
df_regs_full = df_all[(df_all.region > 0.0) & (np.isnan(df_all.lo_mu))]

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# make dfs by mu
plage = df_regs[df_regs.region == 5.0]
network = df_regs[df_regs.region == 4.0]
quiet_sun = df_regs[df_regs.region == 3.0]
penumbrae = df_regs[df_regs.region == 2.0]
umbrae = df_regs[df_regs.region == 1.0]

plage = mask_all_zero_rows(plage)
network = mask_all_zero_rows(network)
quiet_sun = mask_all_zero_rows(quiet_sun)
penumbrae = mask_all_zero_rows(penumbrae)
umbrae = mask_all_zero_rows(umbrae)

# plot the distributions
plt.hist(plage.mag_unsigned, color="tab:purple", bins="auto", histtype="step", density=True)
plt.hist(network.mag_unsigned, color="tab:pink", bins="auto", histtype="step", density=True)
# plt.hist(quiet_sun.mag_unsigned, color="tab:blue", bins="auto", histtype="step", density=True)
plt.hist(penumbrae.mag_unsigned, color="tab:orange", bins="auto", histtype="step", density=True)
plt.hist(umbrae.mag_unsigned, color="tab:green", bins="auto", histtype="step", density=True)

pdb.set_trace()
