import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

# use style
plt.style.use("my.mplstyle"); plt.ioff()

def calc_region_stats(region_df):
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
        reg_avg[i] = np.mean(region_df.v_hat[idx])
        reg_std[i] = np.std(region_df.v_hat[idx])
        reg_err[i] = reg_avg[i]/np.sqrt(len(region_df.v_hat[idx]))
    return reg_avg, reg_std, reg_err

# data files
# datdir = "/Users/michael/Desktop/sdo-pypline/data/"
datdir = "/Users/michael/Desktop/sdo_output/"
df_full = pd.read_csv(datdir + "rv_full_disk.csv")
df_regs = pd.read_csv(datdir + "rv_regions.csv")
df_mu = pd.read_csv(datdir + "rv_mu.csv")

df_full.sort_values("mjd", inplace=True)
df_regs.sort_values("mjd", inplace=True)
df_mu.sort_values("mjd", inplace=True)

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

# get stats
whole_sun_avg, whole_sun_std, whole_sun_err = calc_region_stats(df_mu)
umbrae_avg, umbrae_std, umbrae_err = calc_region_stats(umbrae)
penumbrae_avg, penumbrae_std, penumbrae_err = calc_region_stats(penumbrae)
quiet_sun_avg, quiet_sun_std, quiet_sun_err = calc_region_stats(quiet_sun)
network_avg, network_std, network_err = calc_region_stats(network)
plage_avg, plage_std, plage_err = calc_region_stats(plage)

fig = plt.figure()
ax1 = fig.add_subplot()

ax1.errorbar(mu_bin, whole_sun_avg, yerr=whole_sun_err, fmt=".", color="black", label="All regions")
ax1.fill_between(mu_bin, whole_sun_avg - whole_sun_std, whole_sun_avg + whole_sun_std, color="black", alpha=0.5)

ax1.errorbar(mu_bin, quiet_sun_avg, yerr=quiet_sun_err, fmt=".", color="tab:blue", label="Quiet Sun")
ax1.fill_between(mu_bin, quiet_sun_avg - quiet_sun_std, quiet_sun_avg + quiet_sun_std, color="tab:blue", alpha=0.5)

ax1.errorbar(mu_bin, penumbrae_avg, yerr=penumbrae_err, fmt=".", color="tab:orange", label="Penumbrae")
ax1.fill_between(mu_bin, penumbrae_avg - penumbrae_std, penumbrae_avg + penumbrae_std, color="tab:orange", alpha=0.5)

ax1.errorbar(mu_bin, umbrae_avg, yerr=umbrae_err, fmt=".", color="tab:green", label="Umbrae")
ax1.fill_between(mu_bin, umbrae_avg - umbrae_std, umbrae_avg + umbrae_std, color="tab:green", alpha=0.5)

ax1.errorbar(mu_bin, plage_avg, yerr=plage_err, fmt=".", color="tab:purple", label="Plage")
ax1.fill_between(mu_bin, plage_avg - plage_std, plage_avg + plage_std, color="tab:purple", alpha=0.5)

ax1.errorbar(mu_bin, network_avg, yerr=network_err, fmt=".", color="tab:pink", label="Network")
ax1.fill_between(mu_bin, network_avg - network_std, network_avg + network_std, color="tab:pink", alpha=0.5)

ax1.set_xticks(np.arange(0.1, 1.1, 0.1))
ax1.invert_xaxis()
ax1.set_xlabel(r"$\mu$")
ax1.set_ylabel(r"${\rm Velocity\ (m/s)}$")
ax1.legend(loc="upper left", ncol=3, fontsize=11)
plt.savefig(datdir + "mu_dist_regions.pdf")
plt.clf(); plt.close()
