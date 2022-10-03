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

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# read in the data and sort by mjd
df_all = pd.read_csv(datadir + "velocities.csv")
df_all.sort_values("mjd", inplace=True)

df_stats = pd.read_csv(datadir + "disk_stats.csv")
df_stats.sort_values("mjd", inplace=True)

# get full disk velocities
df_full = df_all[np.isnan(df_all.lo_mu)]

# get velocites for full mu annuli
df_mu = df_all[(df_all.region == 0.0) & (~np.isnan(df_all.lo_mu))]

# get velocities for regions in mu
df_regs = df_all[(df_all.region > 0.0) & (~np.isnan(df_all.lo_mu))]

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.mjd, df_full.v_hat, s=2, label="v_hat")
ax1.scatter(df_full.mjd, df_full.v_phot, s=2, label="v_phot")
ax1.scatter(df_full.mjd, df_full.v_conv, s=2, label="v_conv")
ax1.scatter(df_full.mjd, df_full.v_quiet, s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.legend(ncol=2, fontsize=9)
fig.savefig(plotdir + "full_disk_series.pdf")
plt.clf(); plt.close()

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
idx = (df_regs.region == 1) & (df_regs.v_hat != 0.0)
# ax1.scatter(df_regs.mjd[idx], df_regs.v_hat[idx], s=2, label="v_hat")
ax1.scatter(df_regs.mjd[idx], df_regs.v_phot[idx], s=2, label="v_phot")
ax1.scatter(df_regs.mjd[idx], df_regs.v_conv[idx], s=2, label="v_conv")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_quiet[idx], s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.set_title("Velocity in Spots")
ax1.legend()
fig.savefig(plotdir + "spot_series.pdf")
plt.clf(); plt.close()

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
idx = (df_regs.region == 5) & (df_regs.v_hat != 0.0)
# ax1.scatter(df_regs.mjd[idx], df_regs.v_hat[idx], s=2, label="v_hat")
ax1.scatter(df_regs.mjd[idx], df_regs.v_phot[idx], s=2, label="v_phot")
ax1.scatter(df_regs.mjd[idx], df_regs.v_conv[idx], s=2, label="v_conv")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_quiet[idx], s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.set_title("Velocity in Plage")
ax1.legend()
fig.savefig(plotdir + "plage_series.pdf")
plt.clf(); plt.close()

# time series of other stuff
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_stats.mjd, df_stats.ffactor, s=2, label="Mag. Filling Factor")
ax1.scatter(df_stats.mjd, df_stats.plage_frac, s=2, label="Plage")
ax1.scatter(df_stats.mjd, df_stats.network_frac, s=2, label="Network")
ax1.scatter(df_stats.mjd, df_stats.pen_frac, s=2, label="Penmbrae")
ax1.scatter(df_stats.mjd, df_stats.umb_frac, s=2, label="Umbrae")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_ylim(0.0, 0.12)
ax1.legend(fontsize=9)
fig.savefig(plotdir + "frac_series.pdf")
plt.clf(); plt.close()

# correlation
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_stats.ffactor, df_stats.plage_frac, s=2, label="Plage")
ax1.scatter(df_stats.ffactor, df_stats.plage_frac + df_stats.network_frac, s=2, label="Plage + Network")
ax1.set_xlabel(r"${\rm Magnetic\ Filling\ Factor}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_xlim(0.05, 0.125)
ax1.set_ylim(0.00, 0.125)
ax1.legend(fontsize=9)
fig.savefig(plotdir + "frac_corr.pdf")
plt.clf(); plt.close()
