import numpy as np
import matplotlib.pyplot as plt
import gc, os, re, pdb, csv, glob, time, argparse
from astropy.time import Time
from os.path import exists, split, isdir, getsize

# bring functions into scope
from .paths import root
from .sdo_io import *
from .sdo_vels import *
from .sdo_image import *

# multiprocessing imports
from multiprocessing import get_context
import multiprocessing as mp

def process_data_set_parallel(con_file, mag_file, dop_file,
                              aia_file, mu_thresh, n_rings):
    process_data_set(con_file, mag_file, dop_file, aia_file,
                     mu_thresh=mu_thresh, n_rings=n_rings,
                     suffix=str(mp.current_process().pid))
    return None


def process_data_set(con_file, mag_file, dop_file, aia_file,
                     mu_thresh=0.1, n_rings=10, suffix=None):
    # figure out data directories
    datadir = str(root / "data") + "/"
    if not isdir(datadir):
        os.mkdir(datadir)

    # name output files
    if suffix is None:
        fname1 = datadir + "rv_full_disk.csv"
        fname2 = datadir + "rv_mu.csv"
        fname3 = datadir + "rv_regions.csv"
        fname4 = datadir + "aia_ld_params.csv"
        fname5 = datadir + "hmi_ld_params.csv"
        fname6 = datadir + "con_thresh.csv"
        fname7 = datadir + "mag_stats.csv"
    else:
        # make tmp directory
        tmpdir = datadir + "tmp/"
        if not isdir(tmpdir):
            os.mkdir(tmpdir)

        # filenames
        fname1 = datadir + "tmp/rv_full_disk_" + suffix + ".csv"
        fname2 = datadir + "tmp/rv_mu_" + suffix + ".csv"
        fname3 = datadir + "tmp/rv_regions_" + suffix + ".csv"
        fname4 = datadir + "tmp/aia_ld_params_" + suffix + ".csv"
        fname5 = datadir + "tmp/hmi_ld_params_" + suffix + ".csv"
        fname6 = datadir + "tmp/con_thresh_" + suffix + ".csv"
        fname7 = datadir + "tmp/mag_stats_" + suffix + ".csv"

        # check if the files exist, create otherwise
        if not exists(fname1):
            create_file(fname1)
        if not exists(fname2):
            create_file(fname2)
        if not exists(fname3):
            create_file(fname3)
        if not exists(fname4):
            create_file(fname4)
        if not exists(fname5):
            create_file(fname5)
        if not exists(fname6):
            create_file(fname6)
        if not exists(fname7):
            create_file(fname7)

    # make SDOImage instances
    try:
        con = SDOImage(con_file)
        mag = SDOImage(mag_file)
        dop = SDOImage(dop_file)
        aia = SDOImage(aia_file)
    except OSError:
        print("\t >>> Invalid file, skipping " + get_date(con_file).isoformat(), flush=True)
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
        print("\t >>> Limb darkening fit failed, skipping " + iso, flush=True)
        return None

    # write the limb darkening parameters to disk
    write_results_to_file(fname4, mjd, *aia.ld_coeffs)
    write_results_to_file(fname5, mjd, *con.ld_coeffs)

    # set values to nan for mu less than mu_thresh
    con.mask_low_mu(mu_thresh)
    dop.mask_low_mu(mu_thresh)
    mag.mask_low_mu(mu_thresh)
    aia.mask_low_mu(mu_thresh)

    # identify regions for thresholding
    try:
        mask = SunMask(con, mag, dop, aia)
    except TypeError:
        println("\t >>> Region identification failed, skipping " + iso, flush=True)

    # write thresholds used for masking to disk
    write_results_to_file(fname6, mjd, mask.con_thresh, mask.aia_thresh)

    # compute velocities and write to disk
    vels = calc_velocities(con, mag, dop, aia, mask, region=None, hi_mu=None, lo_mu=None)
    write_results_to_file(fname1, mjd, mask.ff, mask.Bobs, mask.pen_frac,
                          mask.umb_frac, mask.quiet_frac, mask.network_frac,
                          mask.plage_frac, *vels)

    # grid to loop over mu annuli + results containers
    mu_grid = np.linspace(mu_thresh, 1.0, n_rings)
    results_mu = []
    results_reg = []
    results_mag = []

    # calculate disk-integrated unsigned magnetic field
    mag_stats = calc_mag_stats(mag, mask, region=None, lo_mu=None, hi_mu=None)
    results_mag.append([mjd, 0, np.nan, np.nan, *mag_stats])

    # loop over the mu annuli
    for j in range(n_rings-1):
        # mu values for annuli
        lo_mu=mu_grid[j]
        hi_mu=mu_grid[j+1]

        # compute velocity within mu
        vels_mu = calc_velocities(con, mag, dop, aia, mask, region=None, hi_mu=hi_mu, lo_mu=lo_mu)
        results_mu.append([mjd, 0, lo_mu, hi_mu, *vels_mu])

        # calculate disk-integrated unsigned magnetic field
        mag_stats = calc_mag_stats(mag, mask, region=None, hi_mu=hi_mu, lo_mu=lo_mu)
        results_mag.append([mjd, 0, lo_mu, hi_mu, *mag_stats])

        # loop over unique region identifiers
        for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
            # compute velocity components in each mu annulus by region
            vels_reg = calc_velocities(con, mag, dop, aia, mask, region=k, hi_mu=hi_mu, lo_mu=lo_mu)
            results_reg.append([mjd, k, lo_mu, hi_mu, *vels_reg])

            # compute magnetic field strength within each region
            mag_stats = calc_mag_stats(mag, mask, region=k, hi_mu=hi_mu, lo_mu=lo_mu)
            results_mag.append([mjd, k, lo_mu, hi_mu, *mag_stats])

    # write to disk
    write_results_to_file(fname2, results_mu)
    write_results_to_file(fname3, results_reg)
    write_results_to_file(fname7, results_mag)

    # do some memory cleanup
    del con
    del mag
    del dop
    del aia
    del mask
    del vels
    del mu_grid
    del results_mu
    del results_reg
    del results_mag
    gc.collect()

    # report success and return
    print("\t >>> Epoch %s run successfully" % get_date(con_file).isoformat(), flush=True)
    return None
