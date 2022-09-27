import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from sdo_pypline.paths import root
from sdo_pypline.sdo_plot import *
from sdo_pypline.sdo_download import download_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def main():
    # see if data are already downloaded
    files = glob.glob(datadir + "fits/*.fits")
    con = datadir + "fits/hmi_ic_45s_2014_01_07_00_01_30_tai_continuum.fits"
    mag = datadir + "fits/hmi_m_45s_2014_01_07_00_01_30_tai_magnetogram.fits"
    dop = datadir + "fits/hmi_v_45s_2014_01_07_00_01_30_tai_dopplergram.fits"
    aia = datadir + "fits/aia_lev1_1700a_2014_01_07t00_00_30_71z_image_lev1.fits"

    if any(map(lambda x: x not in files, (con, mag, dop, aia))):
        # download the data to plot
        print(">>> Downloading data")
        start = "2014/01/07"
        end = "2014/01/07"
        sample = 24
        con, mag, dop, aia = download_data(outdir=datadir + "fits/", start=start, end=end, sample=sample)

    # preprocess the data and plot it
    print(">>> Processing and plotting data...")
    con = SDOImage(con)
    mag = SDOImage(mag)
    dop = SDOImage(dop)
    aia = SDOImage(aia)

    # calculate geometries
    con.calc_geometry()
    mag.inherit_geometry(con)
    dop.inherit_geometry(con)
    aia.calc_geometry()

    # get MJD for observations
    iso = Time(con.date_obs).iso
    mjd = Time(con.date_obs).mjd

    # interpolate aia image onto hmi image scale
    aia.rescale_to_hmi(con)

    # correct magnetogram for foreshortening
    mag.correct_magnetogram()

    # calculate differential rotation & observer velocity
    dop.calc_vrot_vobs()

    # calculate limb darkening/brightening in continuum map and filtergram
    con.calc_limb_darkening()
    aia.calc_limb_darkening()

    # set values to nan for mu less than mu_thresh
    mu_thresh = 0.1
    con.mask_low_mu(mu_thresh)
    dop.mask_low_mu(mu_thresh)
    mag.mask_low_mu(mu_thresh)
    aia.mask_low_mu(mu_thresh)

    # identify regions for thresholding
    mask = SunMask(con, mag, dop, aia)

    # get intensities to plot
    aia_flat = aia.iflat[~np.isnan(aia.iflat)]
    con_flat = con.iflat[~np.isnan(con.iflat)]
    mag_img = mag.image[~np.isnan(mag.image)]
    w_active = mask.w_active[~np.isnan(mag.image)]
    w_quiet = mask.w_quiet[~np.isnan(mag.image)]

    # length assertion
    assert len(aia_flat) == len(con_flat) == len(mag_img)

    # plot I_flat vs mag for AIA
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axhline(mask.aia_thresh, ls="--", c="k")
    ax1.axvline(24.0, ls="--", c="k")
    ax1.scatter(np.abs(mag_img), aia_flat, s=1, alpha=0.75, color="black", rasterized=True)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\left| B_{r,ij} \right| \ {\rm G}$")
    ax1.set_ylabel(r"${\rm AIA}\ 1700{\rm \AA}\ I_{{\rm flat}, ij}$")
    fig.savefig(plotdir + "aia_flat_vs_mag.pdf", dpi=150)
    plt.clf(); plt.close()

    # plot I_flat vs mag for HMI
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axhline(mask.con_thresh, ls="--", c="k")
    ax1.axvline(24.0, ls="--", c="k")
    ax1.scatter(np.abs(mag_img), con_flat, s=1, alpha=0.75, color="black", rasterized=True)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\left| B_{ij} \right| \ {\rm G}$")
    ax1.set_ylabel(r"${\rm HMI}\ I_{{\rm flat}, ij}$")
    fig.savefig(plotdir + "hmi_flat_vs_mag.pdf", dpi=150)
    plt.clf(); plt.close()

    # plot the distribution of aia intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(mask.aia_thresh, ls="--", c="k")
    _, bins, patches = ax1.hist(aia_flat, cumulative=False, bins="auto", color="black", histtype="step", density=True)
    _, bins, patches = ax1.hist(aia_flat[w_active], cumulative=False, bins="auto", color="tab:blue", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$1700\ {\rm \AA}\ I_{{\rm flat}, ij}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12)
    fig.savefig(plotdir + "aia_flat_dist.pdf")
    plt.clf(); plt.close()

    # plot the distribution of hmi intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(mask.con_thresh, ls="--", c="k")
    _, bins, patches = ax1.hist(con_flat, cumulative=False, bins="auto", color="black", histtype="step", density=True)
    _, bins, patches = ax1.hist(con_flat[w_active], cumulative=False, bins="auto", color="tab:blue", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xlim(3.5e4, ax1.get_xlim()[1])
    ax1.set_xscale("log")
    ax1.set_xlabel(r"${\rm HMI}\ I_{{\rm flat}, ij}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12)
    fig.savefig(plotdir + "hmi_flat_dist.pdf")
    plt.clf(); plt.close()

    # plot the distribution of dark hmi intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(0.6 * mask.con_thresh, ls="dotted", c="k")
    ax1.axvline(mask.con_thresh, ls="--", c="k")
    _, bins, patches = ax1.hist(con_flat[con_flat < mask.con_thresh], cumulative=False, bins="auto", color="black", histtype="step", density=True)
    ax1.set_xscale("log")
    ax1.set_xlabel(r"${\rm HMI}\ I_{{\rm flat}, ij}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    fig.savefig(plotdir + "hmi_flat_dark_dist.pdf")
    plt.clf(); plt.close()

    # plot the distribution of magnetic fields
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.axvline(np.percentile(np.abs(mag_img[con_flat < mask.con_thresh]), 90), ls="dotted", c="k")
    _, bins, patches = ax1.hist(np.abs(mag_img[con_flat < mask.con_thresh]), cumulative=False, bins="auto", color="black", histtype="step", density=True)
    _, bins, patches = ax1.hist(np.abs(mag_img[(con_flat < mask.con_thresh) & (w_active)]), cumulative=False, bins="auto", color="tab:blue", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    _, bins, patches = ax1.hist(np.abs(mag_img[(con_flat < 0.6 * mask.con_thresh) & (w_active)]), cumulative=False, bins="auto", color="tab:orange", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$\left| B_{r,ij} \right| \ {\rm G}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12)
    fig.savefig(plotdir + "mag_dark_dist.pdf")
    plt.clf(); plt.close()

    pdb.set_trace()


    print("derp")
    print("derp")
    print("derp")
    print("derp")
    print("derp")
    print("derp")

    # plot them
    return None

if __name__ == "__main__":
    main()


