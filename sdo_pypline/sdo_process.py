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

def reduce_sdo_images(con_file, mag_file, dop_file, aia_file, mu_thresh=0.1):
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

    # get time of observations
    iso = Time(con.date_obs).iso

    # interpolate aia image onto hmi image scale
    aia.rescale_to_hmi(con)

    # correct magnetogram for foreshortening
    mag.correct_magnetogram()

    # calculate differential rot., meridional circ., obs. vel, grav. redshift, cbs
    dop.correct_dopplergram()

    # calculate limb darkening/brightening in continuum map and filtergram
    try:
        con.calc_limb_darkening()
        aia.calc_limb_darkening()
    except:
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
        mask.mask_low_mu(mu_thresh)
    except:
        print("\t >>> Region identification failed, skipping " + iso, flush=True)
        return None

    return con, mag, dop, aia, mask


def process_data_set_parallel(con_file, mag_file, dop_file, aia_file, mu_thresh, n_rings):
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
        fname1 = datadir + "intensities.csv"
        fname2 = datadir + "pixel_stats.csv"
        fname3 = datadir + "light_stats.csv"
        fname4 = datadir + "velocities.csv"
        fname5 = datadir + "mag_stats.csv"
    else:
        # make tmp directory
        tmpdir = datadir + "tmp/"
        if not isdir(tmpdir):
            os.mkdir(tmpdir)

        # filenames
        fname1 = tmpdir + "intensities_" + suffix + ".csv"
        fname2 = tmpdir + "pixel_stats_" + suffix + ".csv"
        fname3 = tmpdir + "light_stats_" + suffix + ".csv"
        fname4 = tmpdir + "velocities_" + suffix + ".csv"
        fname5 = tmpdir + "mag_stats_" + suffix + ".csv"

        # check if the files exist, create otherwise
        for file in (fname1, fname2, fname3, fname4, fname5):
            if not exists(file):
                create_file(file)

    # reduce the data set
    try:
        con, mag, dop, aia, mask = reduce_sdo_images(con_file, mag_file,
                                                     dop_file, aia_file,
                                                      mu_thresh=mu_thresh)
    except:
        return None

    # get the MJD of the obs
    mjd = Time(con.date_obs).mjd

    # write the limb darkening parameters and pixel fractions to disk
    write_results_to_file(fname1, mjd, mask.aia_thresh, *aia.ld_coeffs,
                          mask.con_thresh1, mask.con_thresh2, *con.ld_coeffs)

    # create arrays to hold velocity magnetic fiel, and pixel fraction results
    results_vel = []
    results_mag = []
    results_pixel = []
    results_light = []

    # append pixel fractions
    all_pixels = np.nansum(con.mu >= mu_thresh)
    results_pixel.append([mjd, np.nan, np.nan, mask.ff, mask.umb_frac,
                         mask.pen_frac, mask.blu_pen_frac, mask.red_pen_frac,
                         mask.quiet_frac, mask.network_frac, mask.plage_frac])

    # calculate light fractions
    all_light = np.nansum(con.image * (con.mu >= mu_thresh))
    umb_light = np.nansum(mask.is_umbra() * con.image)/all_light
    pen_light = np.nansum(mask.is_penumbra() * con.image)/all_light
    blu_pen_light = np.nansum(mask.is_blue_penumbra() * con.image)/all_light
    red_pen_light = np.nansum(mask.is_red_penumbra() * con.image)/all_light
    quiet_light = np.nansum(mask.is_quiet() * con.image)/all_light
    network_light = np.nansum(mask.is_network() * con.image)/all_light
    plage_light = np.nansum(mask.is_plage() * con.image)/all_light

    # append light fractions
    results_pixel.append([mjd, np.nan, np.nan, umb_light, pen_light,
                         blu_pen_light, red_pen_light, quiet_light,
                         network_light, plage_light])


    # calculate disk-integrated velocities
    vels = calc_velocities(con, mag, dop, aia, mask)
    results_vel.append([mjd, 0, np.nan, np.nan, *vels])

    # get v_quiet
    v_quiet = vels[2]

    # calculate disk-integrated unsigned magnetic field
    mags = calc_mag_stats(con, mag, mask)
    results_mag.append([mjd, 0, np.nan, np.nan, *mags])

    # allocate for region mask
    region_mask = np.ones(np.shape(mask.regions)).astype(int)

    # loop over unique region identifiers
    for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
        # compute the region mask
        region_mask[:] = calc_region_mask(mask, region=k, hi_mu=None, lo_mu=None)

        # compute velocity components in each mu annulus by region
        vels = calc_velocities(con, mag, dop, aia, mask, region_mask=region_mask, v_quiet=v_quiet)
        results_vel.append([mjd, k, np.nan, np.nan, *vels])

        # compute magnetic field strength within each region
        mags = calc_mag_stats(con, mag, mask, region_mask=region_mask)
        results_mag.append([mjd, k, np.nan, np.nan, *mags])

    # loop over the mu annuli
    mu_grid = np.linspace(mu_thresh, 1.0, n_rings)
    for j in range(n_rings-1):
        # mu values for annuli
        lo_mu=mu_grid[j]
        hi_mu=mu_grid[j+1]

        # compute the region mask
        region_mask[:] = calc_region_mask(mask, region=None, hi_mu=hi_mu, lo_mu=lo_mu)

        # compute velocity within mu
        vels = calc_velocities(con, mag, dop, aia, mask, region_mask=region_mask)
        results_vel.append([mjd, 0, lo_mu, hi_mu, *vels])

        # get v_quiet for mu annulus
        v_quiet = vels[2]

        # calculate unsigned magnetic field within mu
        mags = calc_mag_stats(con, mag, mask, region_mask=region_mask)
        results_mag.append([mjd, 0, lo_mu, hi_mu, *mags])

        # get filling factor of annulus annulus
        fracs = [np.nansum(mask.w_active * region_mask)/all_pixels]

        # get fraction
        npix_annulus = np.nansum(region_mask)
        fracs = [np.nansum(mask.w_active * region_mask)/npix_annulus]

        # loop over unique region identifiers
        for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
            # compute the region mask
            region_mask[:] = calc_region_mask(mask, region=k, hi_mu=hi_mu, lo_mu=lo_mu)

            # append stats
            fracs.append(np.nansum(region_mask)/npix_annulus)

            # compute velocity components in each mu annulus by region
            vels = calc_velocities(con, mag, dop, aia, mask, region_mask=region_mask, v_quiet=v_quiet)
            results_vel.append([mjd, k, lo_mu, hi_mu, *vels])

            # compute magnetic field strength within each region
            mags = calc_mag_stats(con, mag, mask, region_mask=region_mask)
            results_mag.append([mjd, k, lo_mu, hi_mu, *mags])

        # assemble fractions
        results_pixel.append([mjd, lo_mu, hi_mu, *fracs])

    pdb.set_trace()

    # write to disk
    write_results_to_file(fname2, results_pixel)
    write_results_to_file(fname3, results_light)
    write_results_to_file(fname4, results_vel)
    write_results_to_file(fname5, results_mag)

    # do some memory cleanup
    del con
    del mag
    del dop
    del aia
    del mask
    del vels
    del mags
    del fracs
    del mu_grid
    del region_mask
    del results_vel
    del results_mag
    del results_pixel
    gc.collect()

    # report success and return
    print("\t >>> Epoch %s run successfully" % get_date(con_file).isoformat(), flush=True)
    return None
