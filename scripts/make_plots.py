import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# data files
file_fulldisk = "/Users/michael/Desktop/sdo_output/rv_full_disk.csv"
file_regions = "/Users/michael/Desktop/sdo_output/rv_regions.csv"

# read in
df_full = pd.read_csv(file_fulldisk)
df_regs = pd.read_csv(file_regions)

# sort by date
df_full.sort_values("mjd", inplace=True)
df_regs.sort_values("mjd", inplace=True)

# make time series
plt.scatter(df_full.mjd, df_full.v_hat, label="v_hat")
plt.scatter(df_full.mjd, df_full.v_phot, label="v_phot")
plt.scatter(df_full.mjd, df_full.v_conv, label="v_conv")
plt.scatter(df_full.mjd, df_full.v_quiet, label="v_quiet")
plt.xticks(rotation=45)
plt.xlabel("MJD")
plt.ylabel("Velocity (m/s)")
plt.legend()
plt.savefig("/Users/michael/Desktop/full_disk_series.pdf")
plt.clf(); plt.close()

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# plot by mu
quiet_sun = df_regs[df_regs.region == 3.0]
penumbrae = df_regs[df_regs.region == 1.0]
umbrae = df_regs[df_regs.region == 2.0]

# get mean, std, etc.
quiet_sun_rv_mean = np.zeros(len(lo_mus))
quiet_sun_rv_std = np.zeros(len(lo_mus))
quiet_sun_rv_err = np.zeros(len(lo_mus))
for i in lo_mus:
    idx = quiet_sun.lo_mu == lo_mus[i]
    quiet_sun_rv_mean = np.mean(quiet_sun.v_hat[idx])
    quiet_sun_rv_std = np.std(quiet_sun.v_hat[idx])
    quiet_sun_rv_err = quiet_sun_rv_std/len(quiet_sun.v_hat[idx])

plt.errorbars(mu_bin, quiet_sun_rv_mean, yerr=quiet_sun_rv_err)
plt.xlabel("mu")
plt.ylabel("velocity (m/s)")
plt.savefig("/Users/michael/Desktop/quiet_sun_mu.pdf")


pdb.set_trace()
