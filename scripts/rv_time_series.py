import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# data directories
datdir = "/Users/michael/Desktop/sdo_output/"
df_full = pd.read_csv(datdir + "rv_full_disk.csv")
df_regs = pd.read_csv(datdir + "rv_regions.csv")
df_mu = pd.read_csv(datdir + "rv_mu.csv")

df_full.sort_values("mjd", inplace=True)
df_regs.sort_values("mjd", inplace=True)
df_mu.sort_values("mjd", inplace=True)

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

# # make time series
# fig = plt.figure()
# ax1 = fig.add_subplot()
# idx = (df_regs.region == 1) & (df_regs.v_hat != 0.0)
# # ax1.scatter(df_regs.mjd[idx], df_regs.v_hat[idx], s=2, label="v_hat")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_phot[idx], s=2, label="v_phot")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_conv[idx], s=2, label="v_conv")
# # ax1.scatter(df_regs.mjd[idx], df_regs.v_quiet[idx], s=2, label="v_quiet")
# ax1.set_xlabel(r"${\rm MJD}$")
# ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
# ax1.set_title("Velocity in Spots")
# ax1.legend()
# fig.savefig(datdir + "spot_series.pdf")
# plt.clf(); plt.close()

# # make time series
# fig = plt.figure()
# ax1 = fig.add_subplot()
# idx = (df_regs.region == 5) & (df_regs.v_hat != 0.0)
# # ax1.scatter(df_regs.mjd[idx], df_regs.v_hat[idx], s=2, label="v_hat")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_phot[idx], s=2, label="v_phot")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_conv[idx], s=2, label="v_conv")
# # ax1.scatter(df_regs.mjd[idx], df_regs.v_quiet[idx], s=2, label="v_quiet")
# ax1.set_xlabel(r"${\rm MJD}$")
# ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
# ax1.set_title("Velocity in Plage")
# ax1.legend()
# fig.savefig(datdir + "plage_series.pdf")
# plt.clf(); plt.close()

# time series of other stuff
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.mjd, df_full.ffactor, s=2, label="Mag. Filling Factor")
ax1.scatter(df_full.mjd, df_full.plage_frac, s=2, label="Plage")
ax1.scatter(df_full.mjd, df_full.network_frac, s=2, label="Network")
ax1.scatter(df_full.mjd, df_full.pen_frac, s=2, label="Penmbrae")
ax1.scatter(df_full.mjd, df_full.umb_frac, s=2, label="Umbrae")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_ylim(0.0, 0.12)
ax1.legend(fontsize=9)
fig.savefig(datdir + "frac_series.pdf")
plt.clf(); plt.close()
