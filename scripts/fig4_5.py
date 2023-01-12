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
df_vels = pd.read_csv(datadir + "velocities.csv")
df_vels.sort_values(by=["mjd", "region", "lo_mu"], inplace=True)
df_vels.drop_duplicates()

df_pixels = pd.read_csv(datadir + "pixel_stats.csv")
df_pixels.sort_values(by=["mjd", "lo_mu"], inplace=True)
df_pixels.drop_duplicates()

df_light = pd.read_csv(datadir + "light_stats.csv")
df_light.sort_values(by=["mjd", "lo_mu"], inplace=True)
df_light.drop_duplicates()

df_intensities = pd.read_csv(datadir + "intensities.csv")
df_intensities.sort_values(by=["mjd"], inplace=True)
df_intensities.drop_duplicates()

# get numbers computed for full disk
df_vels_full = df_vels[(np.isnan(df_vels.lo_mu)) & (df_vels.region == 0.0)]
df_pixels_full = df_pixels[np.isnan(df_pixels.lo_mu)]
df_light_full = df_light[np.isnan(df_light.lo_mu)]

# pdb.set_trace()

# fig, ax1 = plt.subplots()
# ax1.scatter(df_vels_full.mjd, df_vels_full.v_conv, s=2, c="k")
# ax1.set_xlabel(r"${\rm MJD}$")
# ax1.set_ylabel(r"$\Delta \hat{v}_{\rm conv}\ {\rm (m/s)}$")
# ax1.set_ylim(0,20)
# fig.savefig("/Users/michael/Desktop/v_conv_time_series.pdf")
# plt.show()


# TODO filtering????
# plt.scatter(df_intensities.mjd, df_intensities.aia_thresh, s=2)
# plt.scatter(df_pixels_full.mjd, df_pixels_full.quiet_frac, s=2)
# plt.scatter(df_pixels_full.mjd, df_pixels_full.plage_frac, s=2)
# plt.scatter(df_pixels_full.mjd, df_pixels_full.network_frac, s=2)
# plt.scatter(df_pixels_full.mjd, df_pixels_full.umb_frac, s=2)
# plt.scatter(df_pixels_full.mjd, df_pixels_full.blu_pen_frac, s=2)
# plt.scatter(df_pixels_full.mjd, df_pixels_full.red_pen_frac, s=2)

# get velocities for regions in mu
df_regs = df_vels[(df_vels.region > 0.0) & (~np.isnan(df_vels.lo_mu))]

# get centers of mu bins
lo_mus = np.unique(df_regs.lo_mu)
hi_mus = np.unique(df_regs.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# make dfs by mu
plage = df_regs[df_regs.region == 6.0]
network = df_regs[df_regs.region == 5.0]
quiet_sun = df_regs[df_regs.region == 4.0]
red_penumbrae = df_regs[df_regs.region == 3.0]
blu_penumbrae = df_regs[df_regs.region == 2.0]
umbrae = df_regs[df_regs.region == 1.0]

# mask rows where all vels are 0.0 (i.e., region isn't present in that annulus)
plage = mask_all_zero_rows(plage)
network = mask_all_zero_rows(network)
quiet_sun = mask_all_zero_rows(quiet_sun)
red_penumbrae = mask_all_zero_rows(red_penumbrae)
blu_penumbrae = mask_all_zero_rows(blu_penumbrae)
umbrae = mask_all_zero_rows(umbrae)


# pdb.set_trace()

# region = quiet_sun[quiet_sun.hi_mu < 0.3]
# plt.scatter(region.mjd, region.v_hat)
# plt.show()


# get stats
def clv_plot(colname, fname=None):
    # whole_avg, whole_std, whole_err = calc_region_stats(df_vels_full, colname=colname)
    umbrae_avg, umbrae_std, umbrae_err = calc_region_stats(umbrae, colname=colname)
    blue_penumbrae_avg, blue_penumbrae_std, blue_penumbrae_err = calc_region_stats(blu_penumbrae, colname=colname)
    red_penumbrae_avg, red_penumbrae_std, red_penumbrae_err = calc_region_stats(red_penumbrae, colname=colname)
    quiet_sun_avg, quiet_sun_std, quiet_sun_err = calc_region_stats(quiet_sun, colname=colname)
    network_avg, network_std, network_err = calc_region_stats(network, colname=colname)
    plage_avg, plage_std, plage_err = calc_region_stats(plage, colname=colname)

    # plot the curves
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(6.4,7.75))
    fig.subplots_adjust(hspace=0.05)

    if colname != "v_conv":
        ax1.errorbar(mu_bin, quiet_sun_avg, yerr=quiet_sun_err, fmt=".", capsize=3, color="tab:blue", label=r"${\rm Quiet\ Sun}$")
        ax1.fill_between(mu_bin, quiet_sun_avg - quiet_sun_std, quiet_sun_avg + quiet_sun_std, color="tab:blue", alpha=0.5)

    # ax1.errorbar(mu_bin, whole_avg, yerr=whole_err, fmt=".", capsize=3, color="k", label=r"${\rm Whole\ Sun}$")
    # ax1.fill_between(mu_bin, whole_avg - whole_std, whole_avg + whole_std, color="k", alpha=0.5)

    ax1.errorbar(mu_bin, plage_avg, yerr=plage_err, fmt=".", capsize=3, color="tab:purple", label=r"${\rm Plage}$")
    ax1.fill_between(mu_bin, plage_avg - plage_std, plage_avg + plage_std, color="tab:purple", alpha=0.5)

    ax1.errorbar(mu_bin, network_avg, yerr=network_err, fmt=".", capsize=3, color="tab:pink", label=r"${\rm Network}$")
    ax1.fill_between(mu_bin, network_avg - network_std, network_avg + network_std, color="tab:pink", alpha=0.5)

    ax2.errorbar(mu_bin, red_penumbrae_avg, yerr=red_penumbrae_err, fmt=".", capsize=3, color="tab:orange", label=r"${\rm Red\ Penumbrae}$")
    ax2.fill_between(mu_bin, red_penumbrae_avg - red_penumbrae_std, red_penumbrae_avg + red_penumbrae_std, color="tab:orange", alpha=0.5)

    ax2.errorbar(mu_bin, blue_penumbrae_avg, yerr=blue_penumbrae_err, fmt=".", capsize=3, color="tab:brown", label=r"${\rm Blue\ Penumbrae}$")
    ax2.fill_between(mu_bin, blue_penumbrae_avg - blue_penumbrae_std, blue_penumbrae_avg + blue_penumbrae_std, color="tab:brown", alpha=0.5)

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
    # if colname == "v_hat":
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
colname = "v_hat"
xlabel = r"$\hat{v}\ {\rm (m/s)}$"
fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx1 = umbrae.lo_mu == mu_samps[i]
    idx2 = red_penumbrae.lo_mu == mu_samps[i]
    idx3 = blu_penumbrae.lo_mu == mu_samps[i]

    # plot this mu
    axs[i].hist(umbrae[colname][idx1], bins="auto", density=True, color="tab:green", histtype="step", label=r"{\rm Umbrae}")
    axs[i].hist(red_penumbrae[colname][idx2], bins="auto", density=True, color="tab:orange", histtype="step", label=r"{\rm Red\ Penumbrae}")
    axs[i].hist(blu_penumbrae[colname][idx3], bins="auto", density=True, color="tab:brown", histtype="step", label=r"{\rm Blue\ Penumbrae}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(xlabel)
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

    # plot this mu
    axs[i].hist(plage[colname][idx1], bins="auto", density=True, color="tab:purple", histtype="step", label=r"{\rm Plage}")
    axs[i].hist(network[colname][idx2], bins="auto", density=True, color="tab:pink", histtype="step", label=r"{\rm Network}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(xlabel)
    axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i] + 0.05)[0:4])
    axs[i].set_xlim(-250,250)
    axs[i].set_box_aspect(1.25)

# set axes labels
axs[0].set_ylabel(r"${\rm Probability\ Density}$")
axs[-1].legend(fontsize=10)
fig.savefig(plotdir + "fig5b.pdf")
plt.clf(); plt.close()

fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx3 = quiet_sun.lo_mu == mu_samps[i]
    axs[i].hist(quiet_sun[colname][idx3], bins="auto", density=True, color="tab:blue", histtype="step", label=r"{\rm Quiet Sun}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(xlabel)
    axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i] + 0.05)[0:4])
    axs[i].set_xlim(-150,150)
    axs[i].set_box_aspect(1.25)
    lbls = axs[i].get_xticklabels()

# set axes labels
axs[0].set_ylabel(r"${\rm Probability\ Density}$")
axs[-1].legend(fontsize=10)
fig.savefig(plotdir + "fig5c.pdf")
plt.clf(); plt.close()


# create figure objects
colname = "v_conv"
xlabel = r"$\Delta \hat{v}_{\rm conv}\ {\rm (m/s)}$"
fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx1 = umbrae.lo_mu == mu_samps[i]
    idx2 = red_penumbrae.lo_mu == mu_samps[i]
    idx3 = blu_penumbrae.lo_mu == mu_samps[i]

    # plot this mu
    axs[i].hist(umbrae[colname][idx1], bins="auto", density=True, color="tab:green", histtype="step", label=r"{\rm Umbrae}")
    axs[i].hist(red_penumbrae[colname][idx2], bins="auto", density=True, color="tab:orange", histtype="step", label=r"{\rm Red\ Penumbrae}")
    axs[i].hist(blu_penumbrae[colname][idx3], bins="auto", density=True, color="tab:brown", histtype="step", label=r"{\rm Blue\ Penumbrae}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(xlabel)
    axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i] + 0.05)[0:4])
    axs[i].set_xlim(-1200, 1200)
    axs[i].set_box_aspect(1.25)

# set axes labels
axs[0].set_ylabel(r"${\rm Probability\ Density}$")
axs[-1].legend(fontsize=10)
fig.savefig(plotdir + "fig6a.pdf")
plt.clf(); plt.close()

fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)

# loop over valus
for i in range(n_mu_samps):
    # do all regs
    idx1 = plage.lo_mu == mu_samps[i]
    idx2 = network.lo_mu == mu_samps[i]

    # plot this mu
    axs[i].hist(plage[colname][idx1], bins="auto", density=True, color="tab:purple", histtype="step", label=r"{\rm Plage}")
    axs[i].hist(network[colname][idx2], bins="auto", density=True, color="tab:pink", histtype="step", label=r"{\rm Network}")

    # plot the full disk
    # axs[i].hist(df_full[colname], bins="auto", density=True, color="k", histtype="step")

    # label stuff
    axs[i].set_xlabel(xlabel)
    axs[i].set_title(r"$\mu =\ $" + str(mu_samps[i] + 0.05)[0:4])
    axs[i].set_xlim(-250,250)
    axs[i].set_box_aspect(1.25)

# set axes labels
axs[0].set_ylabel(r"${\rm Probability\ Density}$")
axs[-1].legend(fontsize=10)
fig.savefig(plotdir + "fig6b.pdf")
plt.clf(); plt.close()

fig, axs = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=n_mu_samps, sharey=True)
fig.subplots_adjust(wspace=0.075)
