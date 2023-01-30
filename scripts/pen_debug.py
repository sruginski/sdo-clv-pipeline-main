import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

from sdo_pypline.paths import root

# use style
plt.style.use("my.mplstyle"); plt.ioff()

def calc_region_stats(region_df, colname="v_hat"):
    # get number elements
    lo_mus = np.unique(region_df.lo_mu[~np.isnan(region_df.lo_mu)])
    nn = len(lo_mus)

    # allocate memory
    reg_avg = np.zeros(nn)
    reg_std = np.zeros(nn)
    reg_err = np.zeros(nn)

    # loop over mu rings
    for i in range(nn):
        # get idx
        idx = region_df.lo_mu == lo_mus[i]

        # calculate the stats
        reg_avg[i] = np.mean(region_df[colname][idx])
        reg_std[i] = np.std(region_df[colname][idx])
        reg_err[i] = reg_avg[i]/np.sqrt(len(region_df[colname][idx]))
    return reg_avg, reg_std, reg_err

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"
procdir = datadir + "processed/"

# read in by region
df_vels_full = pd.read_csv(procdir + "full_disk_vels.csv")
plage = pd.read_csv(procdir + "plage_vels.csv")
network = pd.read_csv(procdir + "network_vels.csv")
quiet_sun = pd.read_csv(procdir + "quiet_sun_vels.csv")
penumbrae = pd.read_csv(procdir + "penumbrae_vels.csv")
red_penumbrae = pd.read_csv(procdir + "red_penumbrae_vels.csv")
blu_penumbrae = pd.read_csv(procdir + "blu_penumbrae_vels.csv")
umbrae = pd.read_csv(procdir + "umbrae_vels.csv")

# read in light
df_light = pd.read_csv(datadir + "light_stats.csv")

# get centers of mu bins
lo_mus = np.unique(plage.lo_mu)
hi_mus = np.unique(plage.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# make plot object and define color rotation
fig, ax1 = plt.subplots(nrows=1, ncols=1)
cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive']

# loop over mu bins
for (j, mu) in enumerate(np.unique(lo_mus)):
    # get penumbrae velocity
    df_temp1 = penumbrae[(penumbrae.lo_mu == mu)]
    df_temp2 = df_light[(df_light.lo_mu == mu) & ((df_light.blu_pen_frac > 0.0) | (df_light.red_pen_frac > 0.0))]

    # throw out dates missing from df_temp1
    common_mjd = pd.merge(df_temp1, df_temp2, how ='inner', on =['mjd']).mjd
    idx = [mjd in common_mjd.values for mjd in df_temp2.mjd.values]
    df_temp2 = df_temp2[idx]

    # get pen light frac
    pen_frac = (df_temp2.blu_pen_frac + df_temp2.red_pen_frac).values# / df_temp2[['umb_frac', 'blu_pen_frac', 'red_pen_frac', 'quiet_frac', 'network_frac', 'plage_frac']].sum(axis=1).values

    # plot it
    # ax1.scatter(df_temp1.v_conv.values, pen_frac, color=cs[j], label=str(mu))

# make the plot pretty
# ax1.set_xlabel("v_conv")
# ax1.set_ylabel("fraction of light")
# ax1.legend()
# fig.savefig("/Users/michael/Desktop/pen_light.pdf")
# plt.show()

