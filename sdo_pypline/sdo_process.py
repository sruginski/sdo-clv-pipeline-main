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

def is_quality_data(sdo_image):
    return sdo_image.quality == 0

def reduce_sdo_images(con_file, mag_file, dop_file, aia_file, mu_thresh=0.1):
    # get the datetime
    iso = get_date(con_file).isoformat()

    # make SDOImage instances
    try:
        con = SDOImage(con_file)
        mag = SDOImage(mag_file)
        dop = SDOImage(dop_file)
        aia = SDOImage(aia_file)
    except OSError:
        print("\t >>> Invalid file, skipping " + iso, flush=True)
        return None

    # check for data quality issue
    if not all(list(map(is_quality_data, [con, mag, dop, aia]))):
        print("\t >>> Data quality issue, skipping " + iso, flush=True)
        return None

    # calculate geometries
    dop.calc_geometry()
    con.inherit_geometry(dop)
    mag.inherit_geometry(dop)

    # interpolate aia image onto hmi image scale and inherit geometry
    aia.rescale_to_hmi(con)

    # calculate limb darkening/brightening in continuum map and filtergram
    try:
        con.calc_limb_darkening()
        aia.calc_limb_darkening()
    except:
        print("\t >>> Limb darkening fit failed, skipping " + iso, flush=True)
        return None

    # correct magnetogram for foreshortening
    mag.correct_magnetogram()

    # calculate differential rot., meridional circ., obs. vel, grav. redshift, cbs
    dop.correct_dopplergram()

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


def process_data_set_parallel(con_file, mag_file, dop_file, aia_file, mu_thresh, n_rings, datadir):
    process_data_set(con_file, mag_file, dop_file, aia_file,
                     mu_thresh=mu_thresh, n_rings=n_rings,
                     suffix=str(mp.current_process().pid), datadir=datadir)
    return None


def process_data_set(con_file, mag_file, dop_file, aia_file,
                     mu_thresh=0.1, n_rings=10, suffix=None, datadir=None):

    # figure out data directories
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
                          mask.blu_pen_frac, mask.red_pen_frac,
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
    results_light.append([mjd, np.nan, np.nan, umb_light, blu_pen_light,
                         red_pen_light, quiet_light, network_light, plage_light])


    # calculate disk-integrated velocities
    vels = calc_velocities(con, mag, dop, aia, mask)
    results_vel.append([mjd, 0, np.nan, np.nan, *vels])

    # calculate disk-integrated unsigned magnetic field
    mags = calc_mag_stats(con, mag, mask)
    results_mag.append([mjd, 0, np.nan, np.nan, *mags])

    # allocate for region mask
    region_mask = np.ones(np.shape(mask.regions)).astype(int)

    # loop over the mu annuli
    mu_grid = np.linspace(mu_thresh, 1.0, n_rings)
    for j in range(n_rings-1):
        # mu values for annuli
        lo_mu=mu_grid[j]
        hi_mu=mu_grid[j+1]

        # compute the region mask
        region_mask[:] = calc_region_mask(mask, region=None, hi_mu=hi_mu, lo_mu=lo_mu)

        # compute quiet-sun velocity in mu annulus
        v_quiet = np.nansum(dop.v_corr * con.image * mask.is_quiet() * region_mask)
        v_quiet /= np.nansum(con.image * mask.is_quiet() * region_mask)

        # get filling factor of annulus
        pixel_fracs = [np.nansum(mask.w_active * region_mask)/all_pixels]
        light_fracs = []

        # loop over unique region identifiers
        for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
            # compute the region mask
            region_mask[:] = calc_region_mask(mask, region=k, hi_mu=hi_mu, lo_mu=lo_mu)

            # calculate and append pixel and light fractions for region
            pixel_fracs.append(np.nansum(region_mask)/all_pixels)
            light_fracs.append(np.nansum(region_mask * con.image)/all_light)

            # compute velocity components in each mu annulus by region
            if k != 4:
                # case where region is quiet sun, return zero for v_quiet
                vels = calc_velocities(con, mag, dop, aia, mask, region_mask=region_mask, v_quiet=v_quiet)
            else:
                # case where region is quiet sun, return nonzero v_quiet
                vels = calc_velocities(con, mag, dop, aia, mask, region_mask=region_mask, v_quiet=None)

            # append the results
            results_vel.append([mjd, k, lo_mu, hi_mu, *vels])

            # compute magnetic field strength within each region
            mags = calc_mag_stats(con, mag, mask, region_mask=region_mask)
            results_mag.append([mjd, k, lo_mu, hi_mu, *mags])

        # assemble fractions
        results_pixel.append([mjd, lo_mu, hi_mu, *pixel_fracs])
        results_light.append([mjd, lo_mu, hi_mu, *light_fracs])

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
    del mu_grid
    del v_quiet
    del pixel_fracs
    del light_fracs
    del region_mask
    del results_vel
    del results_mag
    del results_pixel
    del results_light
    gc.collect()

    # report success and return
    print("\t >>> Epoch %s run successfully" % get_date(con_file).isoformat(), flush=True)
    return None
