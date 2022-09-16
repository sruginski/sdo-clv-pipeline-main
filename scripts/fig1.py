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

    # plot them
    plot_image(con, outdir=plotdir, fname="fig1a.pdf")
    plot_image(mag, outdir=plotdir, fname="fig1b.pdf")
    plot_image(dop, outdir=plotdir, fname="fig1c.pdf")
    plot_image(aia, outdir=plotdir, fname="fig1d.pdf")
    plot_mask(mask, outdir=plotdir, fname="fig2.pdf")
    return None

if __name__ == "__main__":
    main()


