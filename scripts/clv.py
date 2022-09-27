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

def mask_all_zero_rows(df):
    idx = (df.v_hat == 0.0) & (df.v_phot == 0.0) & (df.v_conv == 0.0) & (df.v_quiet == 0.0)
    return df[~idx]

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

# mask rows where all vels are 0.0 (i.e., region isn't present in that area)
plage = mask_all_zero_rows(plage)
network = mask_all_zero_rows(network)
# quiet_sun = mask_all_zero_rows(quiet_sun)
penumbrae = mask_all_zero_rows(penumbrae)
umbrae = mask_all_zero_rows(umbrae)

# get stats
all_regs_avg, all_regs_std, all_regs_err = calc_region_stats(df_mu)
umbrae_avg, umbrae_std, umbrae_err = calc_region_stats(umbrae)
penumbrae_avg, penumbrae_std, penumbrae_err = calc_region_stats(penumbrae)
quiet_sun_avg, quiet_sun_std, quiet_sun_err = calc_region_stats(quiet_sun)
network_avg, network_std, network_err = calc_region_stats(network)
plage_avg, plage_std, plage_err = calc_region_stats(plage)

# plot the curves
fig = plt.figure()
ax1 = fig.add_subplot()

ax1.errorbar(mu_bin, all_regs_avg, yerr=all_regs_err, fmt=".", color="black", label="All regions")
ax1.fill_between(mu_bin, all_regs_avg - all_regs_std, all_regs_avg + all_regs_std, color="black", alpha=0.5)

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

# get distributions of velocities at different mu positions
mu_samps = [0.9, 0.8, 0.4, 0.2]
n_mu_samps = len(mu_samps)

# create figure objects
fig, axs = plt.subplots(figsize=(11, 8.5), nrows=6, ncols=n_mu_samps, tight_layout=True)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx1 = df_mu.lo_mu == mu_samps[i]
    axs[0, i].hist(df_mu.v_hat[idx1], bins="auto", density=True)

    # do the other regions
    idx1 = umbrae.lo_mu == mu_samps[i]
    axs[1, i].hist(umbrae.v_hat[idx1], bins="auto", density=True)

    idx1 = penumbrae.lo_mu == mu_samps[i]
    axs[2, i].hist(penumbrae.v_hat[idx1], bins="auto", density=True)

    idx1 = quiet_sun.lo_mu == mu_samps[i]
    axs[3, i].hist(quiet_sun.v_hat[idx1], bins="auto", density=True)

    idx1 = network.lo_mu == mu_samps[i]
    axs[4, i].hist(network.v_hat[idx1], bins="auto", density=True)

    idx1 = plage.lo_mu == mu_samps[i]
    axs[5, i].hist(plage.v_hat[idx1], bins="auto", density=True)

fig.savefig(datdir + "vel_reg_hists.pdf", bbox_inches="tight")
plt.clf(); plt.close()
