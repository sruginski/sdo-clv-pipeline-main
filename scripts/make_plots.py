import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# data files
# datdir = "/Users/michael/Desktop/sdo-pypline/data/"
datdir = "/Users/michael/Desktop/sdo_output/"
df_full = pd.read_csv(datdir + "rv_full_disk.csv")
df_regs = pd.read_csv(datdir + "rv_regions.csv")
df_mu = pd.read_csv(datdir + "rv_mu.csv")
df_aia_ld = pd.read_csv(datdir + "aia_ld_params.csv")
df_hmi_ld = pd.read_csv(datdir + "hmi_ld_params.csv")
df_con_thresh = pd.read_csv(datdir + "con_thresh.csv")
df_mag_stats = pd.read_csv(datdir + "mag_stats.csv")

# sort by date
df_full.sort_values("mjd", inplace=True)
df_regs.sort_values("mjd", inplace=True)
df_mu.sort_values("mjd", inplace=True)
df_aia_ld.sort_values("mjd", inplace=True)
df_hmi_ld.sort_values("mjd", inplace=True)
df_con_thresh.sort_values("mjd", inplace=True)
df_mag_stats.sort_values("mjd", inplace=True)

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
fig.savefig(datdir + "full_disk_series.pdf")
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
fig.savefig(datdir + "spot_series.pdf")
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
fig.savefig(datdir + "plage_series.pdf")
plt.clf(); plt.close()

# time series of other stuff
network_frac = 1.0 - (df_full.quiet_frac + df_full.pen_frac + df_full.umb_frac + df_full.plage_frac)
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.mjd, df_full.ffactor, s=2, label="Mag. Filling")
ax1.scatter(df_full.mjd, df_full.plage_frac, s=2, label="Plage")
# ax1.scatter(df_full.mjd, df_full.plage_frac + network_frac, s=2, label="Plage + Network")
# ax1.scatter(df_full.mjd, 1.0 - df_full.quiet_frac, s=2, label="Plage + Network + Spot")
ax1.scatter(df_full.mjd, df_full.umb_frac + df_full.pen_frac, s=2, label="Spot")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_ylim(-0.025, 0.12)
ax1.legend(fontsize=9)
fig.savefig(datdir + "frac.pdf")
plt.clf(); plt.close()

# correlation
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.ffactor, df_full.plage_frac, s=2, label="Plage")
ax1.scatter(df_full.ffactor, df_full.plage_frac + network_frac, s=2, label="Plage + Network")
ax1.set_xlabel(r"${\rm Magnetic\ Filling\ Factor}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_xlim(0.05, 0.125)
ax1.set_ylim(0.00, 0.125)
ax1.legend(fontsize=9)
fig.savefig(datdir + "frac_corr.pdf")
plt.clf(); plt.close()

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# plot by mu
plage = df_regs[df_regs.region == 5.0]
network = df_regs[df_regs.region == 4.0]
quiet_sun = df_regs[df_regs.region == 3.0]
penumbrae = df_regs[df_regs.region == 2.0]
umbrae = df_regs[df_regs.region == 1.0]

# plot some histograms
fig = plt.figure()
ax1 = fig.add_subplot()
n, bins, patches = ax1.hist(quiet_sun.v_hat, density=True)
ax1.set_xlabel(r"${\rm Velocity\ (m/s)}$")
ax1.set_ylabel(r"${\rm Probability Density}$")
plt.savefig(datdir + "hist1.pdf")
plt.clf(); plt.close()

# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter()
