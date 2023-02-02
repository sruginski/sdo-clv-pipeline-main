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

def main():
    # get the sdo data to plot
    con, mag, dop, aia = download_plot_data()

    # reduce the data
    print(">>> Processing and plotting data...")
    con, mag, dop, aia, mask = reduce_sdo_images(con, mag, dop, aia)

    pdb.set_trace()

    # plot them
    # print("\t>>> Plotting HMI continuum...")
    plot_image(con, outdir=plotdir, fname="fig1a.pdf")

    # print("\t>>> Plotting magnetogram...")
    plot_image(mag, outdir=plotdir, fname="fig1b.pdf")

    # print("\t>>> Plotting dopplergram...")
    plot_image(dop, outdir=plotdir, fname="fig1c.pdf")

    # print("\t>>> Plotting AIA continuum...")
    plot_image(aia, outdir=plotdir, fname="fig1d.pdf")

    # merge the penumbra
    mask.regions[mask.regions >= 3] -= 1

    # get cutout
    position = (2460, 2250) # in pixels -- origin not lower
    size = u.Quantity((375, 375), u.arcsec)
    cutout_cont = Cutout2D(con.iflat, position, size, wcs=con.wcs)
    cutout_mask = Cutout2D(mask.regions, position, size, wcs=con.wcs)

    # make figure object
    print("\t>>> Plotting mask cutout...")
    fig = plt.figure(figsize=(6.4, 4.8))
    ax1 = fig.add_subplot(111, projection=cutout_cont.wcs)
    ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
    ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
    ax1 = plt.gca() # some weird astropy bug with axes, this fixes it

    # plot a zoom-in on the mask overlaid on continuum
    cmap = plt.get_cmap("Greys").copy()
    cmap.set_bad(color="white")
    ax1.imshow(cutout_cont.data,  cmap=cmap, origin="lower", interpolation=None)

    # cmap = colors.ListedColormap(["black", "saddlebrown", "orange", "yellow", "white"])
    cmap = colors.ListedColormap([um_color, pu_color, qs_color, nw_color, pl_color])
    cmap.set_bad(color="white")
    norm = colors.BoundaryNorm([0, 1, 2, 3, 4, 5], ncolors=cmap.N, clip=True)
    ax1.imshow(cutout_mask.data,  cmap=cmap, alpha=0.5, origin="lower", interpolation=None)

    # axes stuff
    ax1.invert_xaxis()
    ax1.invert_yaxis()
    ax1.grid(False)
    fig.savefig(plotdir + "fig2b.pdf", bbox_inches="tight", dpi=500)
    plt.clf(); plt.close()

    # plot the mask
    print("\t>>> Plotting mask...")
    fig = plt.figure(figsize=(6.4, 4.8))
    ax1 = fig.add_subplot(111, projection=mask.wcs)
    img = ax1.imshow(mask.regions - 0.5, cmap=cmap, norm=norm, origin="lower", interpolation=None)
    sp.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay(ax1, grid_spacing=15*u.deg, annotate=True,
                                                             color="k", alpha=0.5, ls="--", lw=0.5)
    limb = ax1.contour(mask.mu >= 0.0, colors="k", linestyles="--", linewidths=0.5, alpha=0.5)
    cutout_mask.plot_on_original(ax=ax1, color="black", ls="-")
    clb = fig.colorbar(img, ticks=[0.5, 1.5, 2.5, 3.5, 4.5])
    clb.ax.set_yticklabels([r"${\rm Umbra}$", r"${\rm Penumbra}$", r"${\rm Quiet\ Sun}$", r"${\rm Network}$", r"${\rm Plage}$"])
    ax1.invert_xaxis()
    ax1.invert_yaxis()
    ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
    ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
    ax1.text(1400, 4000, mask.date_obs, fontsize=10, c="white")
    ax1.grid(False)
    fig.savefig(plotdir + "fig2a.pdf", bbox_inches="tight", dpi=500)
    plt.clf(); plt.close()

    return None

if __name__ == "__main__":
    main()


