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

def mask_zero_v_conv(df):
    idx = (df.v_hat == df.v_quiet)
    return df[~idx]

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# read in the data and sort by mjd
df_all = pd.read_csv(datadir + "velocities.csv")
df_all.sort_values(by=["mjd"], inplace=True)

df_stats = pd.read_csv(datadir + "disk_stats.csv")
df_stats.sort_values(by=["mjd"], inplace=True)

# get full disk velocities
df_full = df_all[(np.isnan(df_all.lo_mu)) & (df_all.region == 0.0)]

# get velocites for full mu annuli
df_mu = df_all[(~np.isnan(df_all.lo_mu)) & (df_all.region == 0.0)]

# get velocities for regions in mu
df_regs = df_all[(df_all.region > 0.0) & (~np.isnan(df_all.lo_mu))]

# get velocities for regions on fulldisk
df_regs_full = df_all[(df_all.region > 0.0) & (np.isnan(df_all.lo_mu))]

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

# same but whole disk by reg
plage_full = df_regs_full[df_regs_full.region == 5.0]
network_full = df_regs_full[df_regs_full.region == 4.0]
quiet_sun_full = df_regs_full[df_regs_full.region == 3.0]
penumbrae_full = df_regs_full[df_regs_full.region == 2.0]
umbrae_full = df_regs_full[df_regs_full.region == 1.0]

# mask rows where all vels are 0.0 (i.e., region isn't present in that annulus)
plage = mask_all_zero_rows(plage)
network = mask_all_zero_rows(network)
quiet_sun = mask_all_zero_rows(quiet_sun)
penumbrae = mask_all_zero_rows(penumbrae)
umbrae = mask_all_zero_rows(umbrae)

plage_full = mask_all_zero_rows(plage_full)
network_full = mask_all_zero_rows(network_full)
quiet_sun_full = mask_all_zero_rows(quiet_sun_full)
penumbrae_full = mask_all_zero_rows(penumbrae_full)
umbrae_full = mask_all_zero_rows(umbrae_full)

# get stats
def clv_plot(colname, fname=None):
    all_regs_avg, all_regs_std, all_regs_err = calc_region_stats(df_mu, colname=colname)
    umbrae_avg, umbrae_std, umbrae_err = calc_region_stats(umbrae, colname=colname)
    penumbrae_avg, penumbrae_std, penumbrae_err = calc_region_stats(penumbrae, colname=colname)
    quiet_sun_avg, quiet_sun_std, quiet_sun_err = calc_region_stats(quiet_sun, colname=colname)
    network_avg, network_std, network_err = calc_region_stats(network, colname=colname)
    plage_avg, plage_std, plage_err = calc_region_stats(plage, colname=colname)

    # plot the curves
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6.4,7.75))
    fig.subplots_adjust(hspace=0.05)

    ax1.errorbar(mu_bin, all_regs_avg, yerr=all_regs_err, fmt=".", capsize=3, color="black", label=r"${\rm All\ regions}$")
    ax2.errorbar(mu_bin, all_regs_avg, yerr=all_regs_err, fmt=".", capsize=3, color="black")
    ax1.fill_between(mu_bin, all_regs_avg - all_regs_std, all_regs_avg + all_regs_std, color="black", alpha=0.5)
    ax2.fill_between(mu_bin, all_regs_avg - all_regs_std, all_regs_avg + all_regs_std, color="black", alpha=0.5)

    if colname != "v_conv":
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

    if colname == "v_hat":
        ylabel = r"$\hat{v}\ {\rm (m/s)}$"
    elif colname == "v_conv":
        ylabel = r"$\Delta \hat{v}_{\rm conv}\ {\rm (m/s)}$"
    else:
        ylabel = colname

    ax1.set_ylabel(ylabel)
    ax2.set_ylabel(ylabel)

    # stuff for the legend
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    if colname == "v_hat":
        ax1.legend(lines, labels, loc="upper center", ncol=3, fontsize=10)

    # save the figure
    plt.savefig(plotdir + fname, bbox_inches="tight")
    plt.clf(); plt.close()
    return None

clv_plot("v_hat", "fig4a.pdf")
clv_plot("v_conv", "fig4b.pdf")

# get distributions of velocities at different mu positions
mu_samps = [0.9, 0.8, 0.4, 0.2]
n_mu_samps = len(mu_samps)

# create figure objects
colname = "v_conv"
fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx1 = umbrae.lo_mu == mu_samps[i]
    idx2 = penumbrae.lo_mu == mu_samps[i]

    # plot this mu
    axs[i].hist(umbrae[colname][idx1], bins="auto", density=True, color="tab:green", histtype="step", label=r"{\rm Umbrae}")
    axs[i].hist(penumbrae[colname][idx2], bins="auto", density=True, color="tab:orange", histtype="step", label=r"{\rm Penumbrae}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(r"$\hat{v}\ {\rm (m/s)}$")
    axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i] + 0.05)[0:4])
    axs[i].set_xlim(-1200, 1200)
    axs[i].set_box_aspect(1.25)

# set axes labels
axs[0].set_ylabel(r"${\rm Probability\ Density}$")
axs[-1].legend(fontsize=10)
fig.savefig(plotdir + "fig5a.pdf")
plt.clf(); plt.close()

fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx1 = plage.lo_mu == mu_samps[i]
    idx2 = network.lo_mu == mu_samps[i]
    idx3 = quiet_sun.lo_mu == mu_samps[i]

    # plot this mu
    axs[i].hist(plage[colname][idx1], bins="auto", density=True, color="tab:purple", histtype="step", label=r"{\rm Plage}")
    axs[i].hist(network[colname][idx2], bins="auto", density=True, color="tab:pink", histtype="step", label=r"{\rm Network}")
    # axs[i].hist(quiet_sun[colname][idx3], bins="auto", density=True, color="tab:blue", histtype="step", label=r"{\rm Quiet Sun}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(r"$\hat{v}\ {\rm (m/s)}$")
    axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i] + 0.05)[0:4])
    axs[i].set_xlim(-210,450)
    axs[i].set_box_aspect(1.25)

# set axes labels
axs[0].set_ylabel(r"${\rm Probability\ Density}$")
axs[-1].legend(fontsize=10)
fig.savefig(plotdir + "fig5b.pdf")
plt.clf(); plt.close()
