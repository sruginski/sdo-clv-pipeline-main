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

def mask_all_zero_rows(df):
    idx = (df.v_hat == 0.0) & (df.v_phot == 0.0) & (df.v_conv == 0.0) & (df.v_quiet == 0.0)
    return df[~idx]

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
colname = "v_conv"
all_regs_avg, all_regs_std, all_regs_err = calc_region_stats(df_mu, colname=colname)
umbrae_avg, umbrae_std, umbrae_err = calc_region_stats(umbrae, colname=colname)
penumbrae_avg, penumbrae_std, penumbrae_err = calc_region_stats(penumbrae, colname=colname)
quiet_sun_avg, quiet_sun_std, quiet_sun_err = calc_region_stats(quiet_sun, colname=colname)
network_avg, network_std, network_err = calc_region_stats(network, colname=colname)
plage_avg, plage_std, plage_err = calc_region_stats(plage, colname=colname)

# plot the curves
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6.4,7.75))

ax1.errorbar(mu_bin, all_regs_avg, yerr=all_regs_err, fmt=".", capsize=3, color="black", label=r"${\rm All\ regions}$")
ax2.errorbar(mu_bin, all_regs_avg, yerr=all_regs_err, fmt=".", capsize=3, color="black")
ax1.fill_between(mu_bin, all_regs_avg - all_regs_std, all_regs_avg + all_regs_std, color="black", alpha=0.5)
ax2.fill_between(mu_bin, all_regs_avg - all_regs_std, all_regs_avg + all_regs_std, color="black", alpha=0.5)

ax1.errorbar(mu_bin, quiet_sun_avg, yerr=quiet_sun_err, fmt=".", capsize=3, color="tab:blue", label=r"${\rm Quiet\ Sun}$")
ax1.fill_between(mu_bin, quiet_sun_avg - quiet_sun_std, quiet_sun_avg + quiet_sun_std, color="tab:blue", alpha=0.5)

ax1.errorbar(mu_bin, plage_avg, yerr=plage_err, fmt=".", capsize=3, color="tab:purple", label=r"${\rm Plage}$")
ax1.fill_between(mu_bin, plage_avg - plage_std, plage_avg + plage_std, color="tab:purple", alpha=0.5)

ax1.errorbar(mu_bin, network_avg, yerr=network_err, fmt=".", capsize=3, color="tab:pink", label=r"${\rm Network}$")
ax1.fill_between(mu_bin, network_avg - network_std, network_avg + network_std, color="tab:pink", alpha=0.5)

ax2.errorbar(mu_bin, penumbrae_avg, yerr=penumbrae_err, fmt=".", capsize=3, color="tab:orange", label=r"${\rm Penumbrae}$")
ax2.fill_between(mu_bin, penumbrae_avg - penumbrae_std, penumbrae_avg + penumbrae_std, color="tab:orange", alpha=0.5)

ax2.errorbar(mu_bin, umbrae_avg, yerr=umbrae_err, fmt=".", capsize=3, color="tab:green", label=r"${\rm Umbrae}$")
ax2.fill_between(mu_bin, umbrae_avg - umbrae_std, umbrae_avg + umbrae_std, color="tab:green", alpha=0.5)

# annotate axes
ax1.set_xticks(np.arange(0.1, 1.1, 0.1))
ax1.invert_xaxis()
ax2.set_xlabel(r"$\mu$")
ax1.set_ylabel(r"$\hat{v}\ {\rm (m/s)}$")
ax2.set_ylabel(r"$\hat{v}\ {\rm (m/s)}$")

# stuff for the legend
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
ax1.legend(lines, labels, loc="upper center", ncol=3, fontsize=11)

# save the figure
plt.savefig(plotdir + "mu_dist_regions.pdf", bbox_inches="tight")
plt.clf(); plt.close()


# get distributions of velocities at different mu positions
mu_samps = [0.9, 0.8, 0.4, 0.2]
n_mu_samps = len(mu_samps)

def clv_dist_plot(df, color, fname, colname="v_hat", xlims=None):
    # create figure objects
    fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
    fig.subplots_adjust(wspace=0.175, hspace=0.0)

    # loop over valus
    xlims0 = []
    xlims1 = []
    for i in range(n_mu_samps):
        # do all regs
        idx1 = df.lo_mu == mu_samps[i]
        axs[i].hist(df[colname][idx1], bins="auto", density=True, color=color, histtype="step")
        axs[i].set_xlabel(r"$\hat{v}\ {\rm (m/s)}$")
        xlims0.append(axs[i].get_xlim()[0])
        xlims1.append(axs[i].get_xlim()[1])
        axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i]))

    if xlims == None:
        xlim0 = np.min(xlims0)
        xlim1 = np.max(xlims1)
    else:
        xlim0 = xlims[0]
        xlim1 = xlims[1]

    for i in range(n_mu_samps):
        axs[i].set_xlim(xlim0, xlim1)
        axs[i].set_box_aspect(1.25)


    # set axes labels
    axs[0].set_ylabel(r"${\rm Probability\ Density}$")

    fig.savefig(plotdir + fname)
    plt.clf(); plt.close()
    return None

colname = "v_hat"
clv_dist_plot(df_mu, "k", "vel_hist_all.pdf", colname=colname, xlims=(-100,250))
clv_dist_plot(plage, "tab:purple", "vel_hist_plage.pdf", colname=colname, xlims=(-210,450))
clv_dist_plot(network, "tab:pink", "vel_hist_network.pdf", colname=colname, xlims=(-210, 450))
clv_dist_plot(quiet_sun, "tab:blue", "vel_hist_quiet.pdf", colname=colname, xlims=(-100, 260))
clv_dist_plot(penumbrae, "tab:orange", "vel_hist_penumbrae.pdf", colname=colname, xlims=(-1200, 1100))
clv_dist_plot(umbrae, "tab:green", "vel_hist_umbrae.pdf", colname=colname, xlims=(-1200, 1100))
