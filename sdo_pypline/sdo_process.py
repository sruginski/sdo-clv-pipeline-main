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
    else:
        # make tmp directory
        tmpdir = datadir + "tmp/"
        if not isdir(tmpdir):
            os.mkdir(tmpdir)

        # filenames
        fname1 = datadir + "tmp/rv_full_disk_" + suffix + ".csv"
        fname2 = datadir + "tmp/rv_mu_" + suffix + ".csv"
        fname3 = datadir + "tmp/rv_regions_" + suffix + ".csv"

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

    # compute velocities and write to disk
    vels = calc_velocities(con, mag, dop, aia, mask, region=None, hi_mu=None, lo_mu=None)
    write_vels_whole_disk(fname1, mjd, mask.ff, mask.Bobs, mask.pen_frac,
                          mask.umb_frac, mask.quiet_frac, mask.plage_frac, vels)

    # loop over mu annuli
    mu_grid = np.linspace(mu_thresh, 1.0, n_rings)
    results_mu = []
    results_reg = []
    for j in range(n_rings-1):
        # mu values for annuli
        lo_mu=mu_grid[j]
        hi_mu=mu_grid[j+1]

        # compute velocity within mu
        vels_mu = calc_velocities(con, mag, dop, aia, mask, region=None, hi_mu=hi_mu, lo_mu=lo_mu)
        results_mu.append((mjd, 0, lo_mu, hi_mu, *vels_mu))

        # loop over unique region identifiers
        for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
            # compute velocity components in each mu annulus by region
            vels_reg = calc_velocities(con, mag, dop, aia, mask, region=k, hi_mu=hi_mu, lo_mu=lo_mu)
            results_reg.append((mjd, k, lo_mu, hi_mu, *vels_reg))

    # write to disk
    write_vels_by_region(fname2, results_mu)
    write_vels_by_region(fname3, results_reg)

    # do some memory cleanup
    del con
    del mag
    del dop
    del aia
    del mask
    del vels
    del results_mu
    del results_reg
    gc.collect()

    # report success and return
    print("\t >>> Epoch %s run successfully" % get_date(con_file).isoformat(), flush=True)
    return None
