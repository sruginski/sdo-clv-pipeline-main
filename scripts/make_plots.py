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
datdir = "/Users/michael/Desktop/sdo-pypline/data/"
df_full = pd.read_csv(datdir + "rv_full_disk.csv")
df_regs = pd.read_csv(datdir + "rv_regions.csv")
df_mu = pd.read_csv(datdir + "rv_mu.csv")

# sort by date
df_full.sort_values("mjd", inplace=True)
df_regs.sort_values("mjd", inplace=True)
df_mu.sort_values("mjd", inplace=True)

# get coverage of year
# xs = np.arange(np.min(df_full.mjd), np.max(df_full.mjd), np.min(np.diff(df_full.mjd)))
# ys = [any(np.isclose(date, df_full.mjd, rtol=1e-5)) for date in xs]
# plt.scatter(xs, ys, s=1)
# plt.show()

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.mjd, df_full.v_hat, s=2, label="v_hat")
ax1.scatter(df_full.mjd, df_full.v_phot, s=2, label="v_phot")
ax1.scatter(df_full.mjd, df_full.v_conv, s=2, label="v_conv")
ax1.scatter(df_full.mjd, df_full.v_quiet, s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.legend(ncol=4)
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
ax1.scatter(df_full.mjd, df_full.plage_frac + network_frac, s=2, label="Plage + Network")
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

ax1.set_xticks(np.arange(0.1, 1.1, 0.1))
ax1.invert_xaxis()
ax1.set_xlabel(r"$\mu$")
ax1.set_ylabel(r"${\rm Velocity\ (m/s)}$")
ax1.legend(loc="upper left", ncol=3, fontsize=11)
plt.savefig(datdir + "mu_dist_regions.pdf")
plt.clf(); plt.close()



# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter()
