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
file_fulldisk = "/Users/michael/Desktop/sdo_output/rv_full_disk.csv"
file_regions = "/Users/michael/Desktop/sdo_output/rv_regions.csv"
file_mu = "/Users/michael/Desktop/sdo_output/rv_mu.csv"

# read in
df_full = pd.read_csv(file_fulldisk)
df_regs = pd.read_csv(file_regions)
df_mu = pd.read_csv(file_mu)

# sort by date
df_full.sort_values("mjd", inplace=True)
df_regs.sort_values("mjd", inplace=True)
df_mu.sort_values("mjd", inplace=True)

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.mjd, df_full.v_hat, label="v_hat")
ax1.scatter(df_full.mjd, df_full.v_phot, label="v_phot")
ax1.scatter(df_full.mjd, df_full.v_conv, label="v_conv")
ax1.scatter(df_full.mjd, df_full.v_quiet, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.legend()
fig.savefig("/Users/michael/Desktop/full_disk_series.pdf")
plt.clf(); plt.close()

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# plot by mu
plage = df_regs[df_regs.region == 4.0]
quiet_sun = df_regs[df_regs.region == 3.0]
penumbrae = df_regs[df_regs.region == 2.0]
umbrae = df_regs[df_regs.region == 1.0]

# get mean, std, etc.
whole_sun_rv_mean = np.zeros(len(lo_mus))
whole_sun_rv_std = np.zeros(len(lo_mus))
whole_sun_rv_err = np.zeros(len(lo_mus))
quiet_sun_rv_mean = np.zeros(len(lo_mus))
quiet_sun_rv_std = np.zeros(len(lo_mus))
quiet_sun_rv_err = np.zeros(len(lo_mus))
penumbrae_rv_mean = np.zeros(len(lo_mus))
penumbrae_rv_std = np.zeros(len(lo_mus))
penumbrae_rv_err = np.zeros(len(lo_mus))
umbrae_rv_mean = np.zeros(len(lo_mus))
umbrae_rv_std = np.zeros(len(lo_mus))
umbrae_rv_err = np.zeros(len(lo_mus))
plage_rv_mean = np.zeros(len(lo_mus))
plage_rv_std = np.zeros(len(lo_mus))
plage_rv_err = np.zeros(len(lo_mus))
for i in range(len(lo_mus)):
    idx = df_mu.lo_mu == lo_mus[i]
    whole_sun_rv_mean[i] = np.mean(df_mu.v_hat[idx])
    whole_sun_rv_std[i] = np.std(df_mu.v_hat[idx])
    whole_sun_rv_err[i] = whole_sun_rv_mean[i]/np.sqrt(len(df_mu.v_hat[idx]))

    idx = quiet_sun.lo_mu == lo_mus[i]
    quiet_sun_rv_mean[i] = np.mean(quiet_sun.v_hat[idx])
    quiet_sun_rv_std[i] = np.std(quiet_sun.v_hat[idx])
    quiet_sun_rv_err[i] = quiet_sun_rv_std[i]/np.sqrt(len(quiet_sun.v_hat[idx]))

    idx = penumbrae.lo_mu == lo_mus[i]
    penumbrae_rv_mean[i] = np.mean(penumbrae.v_hat[idx])
    penumbrae_rv_std[i] = np.std(penumbrae.v_hat[idx])
    penumbrae_rv_err[i] = penumbrae_rv_std[i]/np.sqrt(len(penumbrae.v_hat[idx]))

    idx = umbrae.lo_mu == lo_mus[i]
    umbrae_rv_mean[i] = np.mean(umbrae.v_hat[idx])
    umbrae_rv_std[i] = np.std(umbrae.v_hat[idx])
    umbrae_rv_err[i] = umbrae_rv_std[i]/np.sqrt(len(umbrae.v_hat[idx]))

    idx = plage.lo_mu == lo_mus[i]
    plage_rv_mean[i] = np.mean(plage.v_hat[idx])
    plage_rv_std[i] = np.std(plage.v_hat[idx])
    plage_rv_err[i] = plage_rv_std[i]/np.sqrt(len(plage.v_hat[idx]))


fig = plt.figure()
ax1 = fig.add_subplot()

ax1.errorbar(mu_bin, whole_sun_rv_mean, yerr=whole_sun_rv_err, fmt=".", color="black", label="All regions")
ax1.fill_between(mu_bin, whole_sun_rv_mean - whole_sun_rv_std, whole_sun_rv_mean + whole_sun_rv_std, color="black", alpha=0.5)

ax1.errorbar(mu_bin, quiet_sun_rv_mean, yerr=quiet_sun_rv_err, fmt=".", color="tab:blue", label="Quiet Sun")
ax1.fill_between(mu_bin, quiet_sun_rv_mean - quiet_sun_rv_std, quiet_sun_rv_mean + quiet_sun_rv_std, color="tab:blue", alpha=0.5)

# ax1.errorbar(mu_bin, penumbrae_rv_mean, yerr=penumbrae_rv_err, fmt=".", color="tab:orange", label="Penumbrae")
# ax1.fill_between(mu_bin, penumbrae_rv_mean - penumbrae_rv_std, penumbrae_rv_mean + penumbrae_rv_std, color="tab:orange", alpha=0.5)

# ax1.errorbar(mu_bin, umbrae_rv_mean, yerr=umbrae_rv_err, fmt=".", color="tab:green", label="Umbrae")
# ax1.fill_between(mu_bin, umbrae_rv_mean - umbrae_rv_std, umbrae_rv_mean + umbrae_rv_std, color="tab:green", alpha=0.5)

ax1.errorbar(mu_bin, plage_rv_mean, yerr=plage_rv_err, fmt=".", color="tab:purple", label="Plage")
ax1.fill_between(mu_bin, plage_rv_mean - plage_rv_std, plage_rv_mean + plage_rv_std, color="tab:purple", alpha=0.5)

ax1.invert_xaxis()
ax1.set_xlabel(r"$\mu$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.legend(loc="upper left")
plt.savefig("/Users/michael/Desktop/mu_dist_regions.pdf")
plt.clf(); plt.close()
