import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from astropy.nddata import Cutout2D
from astropy import units as u

from sdo_pypline.paths import root
from sdo_pypline.sdo_plot import *
from sdo_pypline.sdo_process import *
from download_plot_data import download_plot_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

pl_color = "tab:purple" # colors[4]
nw_color = "tab:pink" # colors[6]
qs_color = "orange" # colors[1]
rp_color = "tab:red" # colors[3]
bp_color = "tab:blue" # colors[9]
pu_color = "sienna" # colors[5]
um_color = "tab:gray" # colors[7]

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"
procdir = datadir + "processed/"

# read penumbrae data
penumbrae = pd.read_csv(procdir + "penumbrae.csv")
red_penumbrae = pd.read_csv(procdir + "red_penumbrae.csv")
blu_penumbrae = pd.read_csv(procdir + "blu_penumbrae.csv")

# get the LD coeffs by day
thresholds = pd.read_csv(datadir + "thresholds.csv")
thresholds.sort_values(by=["mjd"], inplace=True)

# get centers of mu bins
lo_mus = np.unique(penumbrae.lo_mu)
hi_mus = np.unique(penumbrae.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# get avg ints
ap_ints = np.zeros(len(penumbrae.avg_int))
rp_ints = np.zeros(len(red_penumbrae.avg_int))
bp_ints = np.zeros(len(blu_penumbrae.avg_int))
ap_ints_flat = np.zeros(len(penumbrae.avg_int))
rp_ints_flat = np.zeros(len(red_penumbrae.avg_int))
bp_ints_flat = np.zeros(len(blu_penumbrae.avg_int))

# loop over mjd and divide out 0th LD parameter
for (i, mjd) in enumerate(thresholds.mjd):
    idx1 = penumbrae.mjd == mjd
    idx2 = red_penumbrae.mjd == mjd
    idx3 = blu_penumbrae.mjd == mjd
    ap_ints[np.where(idx1)] = penumbrae.avg_int[idx1].values/thresholds.a_hmi[thresholds.mjd == mjd].values
    rp_ints[np.where(idx2)] = red_penumbrae.avg_int[idx2].values/thresholds.a_hmi[thresholds.mjd == mjd].values
    bp_ints[np.where(idx3)] = blu_penumbrae.avg_int[idx3].values/thresholds.a_hmi[thresholds.mjd == mjd].values
    ap_ints_flat[np.where(idx1)] = penumbrae.avg_int_flat[idx1].values/thresholds.a_hmi[thresholds.mjd == mjd].values
    rp_ints_flat[np.where(idx2)] = red_penumbrae.avg_int_flat[idx2].values/thresholds.a_hmi[thresholds.mjd == mjd].values
    bp_ints_flat[np.where(idx3)] = blu_penumbrae.avg_int_flat[idx3].values/thresholds.a_hmi[thresholds.mjd == mjd].values

# plot it
fig, ax1 = plt.subplots()
ax1.hist(ap_ints_flat, density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
ax1.hist(rp_ints_flat, density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
ax1.hist(bp_ints_flat, density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")
ax1.set_xlabel(r"{\rm Relative\ Flattened\ Intensity}")
ax1.set_ylabel(r"{\rm Probability\ Density}")
ax1.set_xlim(0.59, 0.91)
ax1.legend(loc="upper left")
fig.savefig(plotdir + "fig7a.pdf")
plt.clf(); plt.close()

# plot by mu
fig, axs = plt.subplots(figsize=(8.75, 7), nrows=3, ncols=3, sharey=True, sharex=True)
fig.subplots_adjust(wspace=0.05)
for (i, mu) in enumerate(lo_mus):
    # get indices for data
    idx1 = np.where(penumbrae.lo_mu == mu)
    idx2 = np.where(red_penumbrae.lo_mu == mu)
    idx3 = np.where(blu_penumbrae.lo_mu == mu)

    # plot it
    axn = fig.axes[i]
    axn.hist(ap_ints[idx1], density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
    axn.hist(rp_ints[idx2], density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
    axn.hist(bp_ints[idx3], density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")

    # set the title
    axn.set_title(r"$\mu = " + str(mu + 0.05)[0:4]+ r"$")

# label the axes
fig.supxlabel(r"{\rm Relative\ Intensity}")
fig.supylabel(r"{\rm Probability\ Density}")

# prepare the legend
handles = [axn.get_legend_handles_labels()[0] for axn in fig.axes][0]
labels = [axn.get_legend_handles_labels()[1] for axn in fig.axes][0]
fig.legend(handles, labels, ncol=3, loc='upper center', handletextpad=0.15, bbox_to_anchor=(0.51, 1.0))
fig.savefig(plotdir + "fig7b.pdf", bbox_inches="tight")
plt.clf(); plt.close()


# plot by mu
fig, axs = plt.subplots(figsize=(8.75, 7), nrows=3, ncols=3, sharey=True, sharex=True)
fig.subplots_adjust(wspace=0.05)
for (i, mu) in enumerate(lo_mus):
    # get indices for data
    idx1 = np.where(penumbrae.lo_mu == mu)
    idx2 = np.where(red_penumbrae.lo_mu == mu)
    idx3 = np.where(blu_penumbrae.lo_mu == mu)

    # plot it
    axn = fig.axes[i]
    axn.hist(ap_ints_flat[idx1], density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
    axn.hist(rp_ints_flat[idx2], density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
    axn.hist(bp_ints_flat[idx3], density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")

    # set the title
    axn.set_title(r"$\mu = " + str(mu + 0.05)[0:4]+ r"$")

# label the axes
fig.supxlabel(r"{\rm Relative\ Flattened\ Intensity}")
fig.supylabel(r"{\rm Probability\ Density}")

# prepare the legend
handles = [axn.get_legend_handles_labels()[0] for axn in fig.axes][0]
labels = [axn.get_legend_handles_labels()[1] for axn in fig.axes][0]
fig.legend(handles, labels, ncol=3, loc='upper center', handletextpad=0.15, bbox_to_anchor=(0.51, 1.0))
fig.savefig(plotdir + "fig7c.pdf", bbox_inches="tight")
plt.clf(); plt.close()
