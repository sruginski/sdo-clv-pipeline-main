import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import re, pdb, csv, glob, time, argparse
from astropy.time import Time
from os.path import exists, split, isdir, getsize

from .paths import root
from .sdo_io import *
from .sdo_plot import *
from .sdo_vels import *
from .sdo_image import *


def process_data_set(con_file, mag_file, dop_file, aia_file, mu_thresh=0.1,
                     n_rings=10, plot=False, vels=True, **kwargs):
    # figure out data directories
    if "datadir" not in kwargs:
        datadir = str(root / "data") + "/"
    else:
        datadir = kwargs["datadir"]

    # name output files
    fname1 = datadir + "rv_full_disk.csv"
    fname2 = datadir + "rv_mu.csv"
    fname3 = datadir + "rv_regions.csv"

    # figure out plot output directors
    if "plotdir" not in kwargs:
        plotdir = str(root / "figures") + "/"
    else:
        plotdir = kwargs["plotdir"]

    # make SDOImage instances
    try:
        con = SDOImage(con_file)
        mag = SDOImage(mag_file)
        dop = SDOImage(dop_file)
        aia = SDOImage(aia_file)
    except OSError:
        print("\t >>> Invalid file, skipping " + get_date(con_file).isoformat())
        return None

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
    try:
        con.calc_limb_darkening()
        aia.calc_limb_darkening()
    except RuntimeError:
        print("\t >>> Limb darkening fit failed, skipping " + iso)
        return None

    # set values to nan for mu less than mu_thresh
    con.mask_low_mu(mu_thresh)
    dop.mask_low_mu(mu_thresh)
    mag.mask_low_mu(mu_thresh)
    aia.mask_low_mu(mu_thresh)

    # identify regions for thresholding
    mask = SunMask(con, mag, dop, aia)

    # plot the data
    if plot:
        plot_image(con, outdir=plotdir, **kwargs)
        plot_image(mag, outdir=plotdir, **kwargs)
        plot_image(dop, outdir=plotdir, **kwargs)
        plot_image(aia, outdir=plotdir, **kwargs)
        plot_mask(mask, outdir=plotdir, **kwargs)

    if vels:
        # compute velocities and write to disk
        vels = calc_velocities(con, mag, dop, aia, mask)
        write_vels(fname1, mjd, mask.ff, mask.Bobs, mask.pen_frac,
                   mask.umb_frac, mask.quiet_frac,
                   mask.plage_frac, vels)

        # loop over mu annuli
        mu_grid = np.linspace(mu_thresh, 1.0, n_rings)
        for j in range(n_rings-1):
            # mu values for annuli
            lo_mu=mu_grid[j]
            hi_mu=mu_grid[j+1]

            # # compute velocity in mu annulus
            # vels_reg = calc_velocities(con, mag, dop, aia, mask,
            #                            lo_mu=lo_mu, hi_mu=hi_mu)

            # # write to disk
            # write_vels_by_region(fname2, mjd, 0, lo_mu, hi_mu, vels_reg)

            # loop over unique region identifiers
            for j in np.unique(mask.regions[~np.isnan(mask.regions)]):
                # compute velocity components in each mu annulus by region
                vels_reg = calc_velocities(con, mag, dop, aia, mask, region=j,
                                           lo_mu=lo_mu, hi_mu=hi_mu)

                # write to disk
                write_vels_by_region(fname3, mjd, j, lo_mu, hi_mu, vels_reg)

    return None
