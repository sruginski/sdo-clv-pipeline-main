import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from sdo_pypline.paths import root
from sdo_pypline.sdo_process import *
from download_plot_data import download_plot_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def main():
    # get the sdo data to plot
    con, mag, dop, aia = download_plot_data()

    # reduce the data
    print(">>> Processing and plotting data...")
    con, mag, dop, aia, mask = reduce_sdo_images(con, mag, dop, aia, mu_thresh=0.1)

    # get intensities to plot
    aia_flat = aia.iflat[~np.isnan(aia.iflat)] / aia.ld_coeffs[0]
    con_flat = con.iflat[~np.isnan(con.iflat)] / con.ld_coeffs[0]
    mag_img = mag.image[~np.isnan(mag.image)]
    w_active = mask.w_active[~np.isnan(mag.image)]
    w_quiet = mask.w_quiet[~np.isnan(mag.image)]
    is_umbra = mask.is_umbra()[~np.isnan(mag.image)]
    is_penumbra = mask.is_penumbra()[~np.isnan(mag.image)]
    is_dark = np.logical_or(mask.is_umbra(), mask.is_penumbra())[~np.isnan(mag.image)]
    is_network = mask.is_network()[~np.isnan(mag.image)]
    is_plage = mask.is_plage()[~np.isnan(mag.image)]
    is_bright = np.logical_or(mask.is_network(), mask.is_plage())[~np.isnan(mag.image)]
    mask_dark = ~is_dark

    # length assertion
    assert len(aia_flat) == len(con_flat) == len(mag_img)

    # get random indices
    idx1 = np.random.choice(len(aia_flat), int(np.floor(len(aia_flat)/50)))
    idx2 = np.random.choice(len(aia_flat[w_active]), int(np.floor(len(aia_flat[w_active])/50)))

    # plot the distribution of intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(mask.con_thresh1/con.ld_coeffs[0], ls="--", c="k", label=r"$0.89\ \hat{I}_{\rm quiet}$")
    ax1.axvspan(0, mask.con_thresh1/con.ld_coeffs[0], fill=False, ec="tab:orange", hatch="/", alpha=0.75)
    ax1.axvspan(0, mask.con_thresh1/con.ld_coeffs[0], fill=False, ec="tab:green", hatch="\\", alpha=0.75)
    x1, bins, patches = ax1.hist(con_flat, bins="fd", color="black", lw=2, histtype="step", density=True, label=r"${\rm Any\ }\left| B_{r}\right|$")
    x2, bins, patches = ax1.hist(con_flat[w_active], bins="fd", color="tab:blue", lw=2, histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlim(6e-1, 1.25)
    ax1.set_xlabel(r"${\rm Normalized\ HMI\ Continuum\ Intensity}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12, loc="upper left")
    ax1.grid(False)
    fig.savefig(plotdir + "fig3a.pdf")
    plt.clf(); plt.close()

    # plot the distribution of dark hmi intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(mask.con_thresh1/con.ld_coeffs[0], ls="--", c="k", label=r"$0.89\ \hat{I}_{\rm quiet}$")
    ax1.axvline(mask.con_thresh2/con.ld_coeffs[0], ls=":", c="k", label=r"$0.45\ \hat{I}_{\rm quiet}$")
    ax1.axvspan(0, mask.con_thresh2/con.ld_coeffs[0], fill=False, ec="tab:green", hatch="/", alpha=0.75)
    ax1.axvspan(mask.con_thresh2/con.ld_coeffs[0], mask.con_thresh1/con.ld_coeffs[0], fill=False, ec="tab:orange", hatch="\\", alpha=0.75)
    x1, bins, patches = ax1.hist(con_flat[con_flat < (mask.con_thresh1/con.ld_coeffs[0])], cumulative=False, bins="auto", lw=2, color="black", histtype="step", density=True, label=r"${\rm Any\ }\left| B_{r}\right|$")
    x2, bins, patches = ax1.hist(con_flat[(con_flat < (mask.con_thresh1/con.ld_coeffs[0])) & w_active], cumulative=False, bins="auto", lw=2, color="tab:blue", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"${\rm Normalized\ HMI\ Continuum\ Intensity}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12, loc="upper left")
    ax1.grid(False)
    fig.savefig(plotdir + "fig3b.pdf")
    plt.clf(); plt.close()

    # plot them
    return None

if __name__ == "__main__":
    main()


