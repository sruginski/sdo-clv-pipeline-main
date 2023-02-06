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
quiet_sun = pd.read_csv(procdir + "quiet_sun.csv")

# get the LD coeffs by day
thresholds = pd.read_csv(datadir + "thresholds.csv")

# get centers of mu bins
lo_mus = np.unique(penumbrae.lo_mu)
hi_mus = np.unique(penumbrae.hi_mu)
mu_bin = (lo_mus + hi_mus) / 2.0

# get avg ints aggregated by mu
ap = penumbrae.avg_int_flat
rp = red_penumbrae.avg_int_flat
bp = blu_penumbrae.avg_int_flat
qs = quiet_sun.avg_int_flat

# plot it
fig, ax1 = plt.subplots()
ax1.hist(ap, density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
ax1.hist(rp, density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
ax1.hist(bp, density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")
ax1.set_xlabel(r"{\rm Flattened\ Intensity}")
ax1.set_ylabel(r"{\rm Probability\ Density}")
ax1.legend(loc="upper left")
fig.savefig(plotdir + "fig7a.pdf")
plt.clf(); plt.close()

# plot by mu
fig, axs = plt.subplots(figsize=(8.75, 7), nrows=3, ncols=3, sharey=True, sharex=True)
# fig, axs = plt.subplots(nrows=3, ncols=3, sharey=True, sharex=True)
fig.subplots_adjust(wspace=0.05)
for (i, mu) in enumerate(lo_mus):
    # get temp data frames
    ap = penumbrae[penumbrae.lo_mu == mu].avg_int_flat
    rp = red_penumbrae[red_penumbrae.lo_mu == mu].avg_int_flat
    bp = blu_penumbrae[blu_penumbrae.lo_mu == mu].avg_int_flat

    # plot it
    axn = fig.axes[i]
    axn.hist(ap, density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
    axn.hist(rp, density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
    axn.hist(bp, density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")

    # set the title
    axn.set_title(r"$\mu = " + str(mu + 0.05)[0:4]+ r"$")

# label the axes
fig.supxlabel(r"{\rm Flattened\ Intensity}")
fig.supylabel(r"{\rm Probability\ Density}")

# prepare the legend
handles = [axn.get_legend_handles_labels()[0] for axn in fig.axes][0]
labels = [axn.get_legend_handles_labels()[1] for axn in fig.axes][0]
fig.legend(handles, labels, ncol=3, loc='upper center', handletextpad=0.15, bbox_to_anchor=(0.51, 1.0))
plt.show()


# # get the sdo data to plot
# con, mag, dop, aia = download_plot_data()

# # reduce the data
# print(">>> Processing and plotting data...")
# con, mag, dop, aia, mask = reduce_sdo_images(con, mag, dop, aia)

# # get masks
# red_pen_mask = mask.is_red_penumbra()
# blu_pen_mask = mask.is_blue_penumbra()
# all_pen_mask = mask.is_penumbra()

# # get intensities
# red_pen_int = con.image[red_pen_mask]
# blu_pen_int = con.image[blu_pen_mask]
# all_pen_int = con.image[all_pen_mask]
# red_pen_int_flat = con.iflat[red_pen_mask]
# blu_pen_int_flat = con.iflat[blu_pen_mask]
# all_pen_int_flat = con.iflat[all_pen_mask]

# # get velocities
# red_pen_vel = dop.v_corr[red_pen_mask]
# blu_pen_vel = dop.v_corr[blu_pen_mask]
# all_pen_vel = dop.v_corr[all_pen_mask]

# # scatter plot of intensities
# fig, ax1 = plt.subplots()
# ax1.scatter(red_pen_int/con.ld_coeffs[0], np.abs(red_pen_vel), marker="v", c=rp_color, alpha=0.5, label=r"${\rm Red\ Penumbrae}$", rasterized=True)
# ax1.scatter(blu_pen_int/con.ld_coeffs[0], np.abs(blu_pen_vel), marker="^", c=bp_color, alpha=0.5, label=r"${\rm Blue\ Penumbrae}$", rasterized=True)
# ax1.set_xlabel(r"${\rm Relative\ Intensity}$")
# ax1.set_ylabel(r"$\left| {\rm Velocity\ } \right| \ {\rm(m\ s}^{-1}{\rm )}$")
# ax1.legend()
# fig.savefig("/Users/michael/Desktop/intensity_velocity.pdf", dpi=250)
# plt.clf(); plt.close()

# fig, ax1 = plt.subplots()
# ax1.scatter(red_pen_int_flat/con.ld_coeffs[0], np.abs(red_pen_vel), marker="v", c=rp_color, alpha=0.5, label=r"${\rm Red\ Penumbrae}$", rasterized=True)
# ax1.scatter(blu_pen_int_flat/con.ld_coeffs[0], np.abs(blu_pen_vel), marker="^", c=bp_color, alpha=0.5, label=r"${\rm Blue\ Penumbrae}$", rasterized=True)
# ax1.set_xlabel(r"${\rm Relative\ Flattened\ Intensity}$")
# ax1.set_ylabel(r"$\left| {\rm Velocity\ } \right| \ {\rm(m\ s}^{-1}{\rm )}$")
# ax1.legend()
# fig.savefig("/Users/michael/Desktop/flattened_intensity_velocity.pdf", dpi=250)
# plt.clf(); plt.close()

# # distributions of intensities
# fig, ax1 = plt.subplots()
# ax1.hist(red_pen_int/con.ld_coeffs[0], density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
# ax1.hist(blu_pen_int/con.ld_coeffs[0], density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")
# ax1.hist(all_pen_int/con.ld_coeffs[0], density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
# ax1.set_xlabel(r"${\rm Relative\ Intensity}$")
# ax1.set_ylabel(r"${\rm Probability\ Density}$")
# ax1.legend(loc="upper left")
# fig.savefig("/Users/michael/Desktop/intensity_distribution.pdf")
# plt.clf(); plt.close()

# fig, ax1 = plt.subplots()
# ax1.hist(red_pen_int_flat/con.ld_coeffs[0], density=True, bins="auto", histtype="step", color=rp_color, label=r"${\rm Red\ Penumbrae}$")
# ax1.hist(blu_pen_int_flat/con.ld_coeffs[0], density=True, bins="auto", histtype="step", color=bp_color, label=r"${\rm Blue\ Penumbrae}$")
# ax1.hist(all_pen_int_flat/con.ld_coeffs[0], density=True, bins="auto", histtype="step", color=pu_color, label=r"${\rm All\ Penumbrae}$")
# ax1.set_xlabel(r"${\rm Relative\ Flattened\ Intensity}$")
# ax1.set_ylabel(r"${\rm Probability\ Density}$")
# ax1.legend(loc="upper left")
# fig.savefig("/Users/michael/Desktop/flattened_intensity_distribution.pdf")
# plt.clf(); plt.close()

# pdb.set_trace()
