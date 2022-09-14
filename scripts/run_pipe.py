import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import os, re, pdb, csv, glob, time, argparse
from astropy.time import Time
from os.path import exists, split, isdir, getsize

# bring functions into scope
from sdo_pypline.paths import root
from sdo_pypline.sdo_io import *
from sdo_pypline.sdo_vels import *
from sdo_pypline.sdo_plot import *
from sdo_pypline.sdo_image import *

# multiprocessing imports
from multiprocessing import get_context
import multiprocessing as mp

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def process_data_set_parallel(con_file, mag_file, dop_file,
                              aia_file, mu_thresh, n_rings):
    process_data_set(con_file, mag_file, dop_file, aia_file,
                     mu_thresh=mu_thresh, n_rings=n_rings,
                     suffix=str(mp.current_process().pid),
                     plot=False, vels=True)
    return None

def process_data_set(con_file, mag_file, dop_file, aia_file, mu_thresh=0.1,
                     n_rings=10, plot=False, vels=True, suffix=None, **kwargs):
    # figure out data directories
    if "datadir" not in kwargs:
        datadir = str(root / "data") + "/"
    else:
        datadir = kwargs["datadir"]

    # name output files
    if suffix is None:
        fname1 = datadir + "rv_full_disk.csv"
        fname2 = datadir + "rv_mu.csv"
        fname3 = datadir + "rv_regions.csv"
    else:
        fname1 = datadir + "rv_full_disk_" + suffix + ".csv"
        fname2 = datadir + "rv_mu_" + suffix + ".csv"
        fname3 = datadir + "rv_regions_" + suffix + ".csv"

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
        vels = calc_velocities(con, mag, dop, aia, mask, None, None, None)
        write_vels_whole_disk(fname1, mjd, mask.ff, mask.Bobs, mask.pen_frac,
                              mask.umb_frac, mask.quiet_frac, mask.plage_frac, vels)

        # loop over mu annuli
        mu_grid = np.linspace(mu_thresh, 1.0, n_rings)
        results = []
        for j in range(n_rings-1):
            # mu values for annuli
            lo_mu=mu_grid[j]
            hi_mu=mu_grid[j+1]

            # loop over unique region identifiers
            for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
                # compute velocity components in each mu annulus by region
                vels_reg = calc_velocities(con, mag, dop, aia, mask, k, hi_mu, lo_mu)
                results.append((mjd, k, lo_mu, hi_mu, *vels_reg))

        # write to disk
        write_vels_by_region(fname3, results)
    print("\t >>> Epoch %s run successfully" % get_date(con_file).isoformat(), flush=True)
    return None

def get_parser_args():
    # initialize argparser
    parser = argparse.ArgumentParser(description="Analyze SDO data")
    parser.add_argument("--clobber", action="store_true", default=False)

    # parse the command line arguments
    args = parser.parse_args()
    clobber = args.clobber
    return clobber

def main():
    # define sdo_data directories
    indir = "/Users/michael/Desktop/sdo_data/"
    if not isdir(indir):
        indir = "/storage/home/mlp95/scratch/sdo_data/"

    # sort out input/output data files
    clobber = get_parser_args()
    con_files, mag_files, dop_files, aia_files = organize_input_output(indir, clobber=clobber)

    # set mu threshold, number of mu rings
    n_rings = 10
    mu_thresh = 0.1
    plot = False

    # get number of cpus
    try:
        from os import sched_getaffinity
        ncpus = len(sched_getaffinity(0))
    except:
        # ncpus = np.min([len(con_files), mp.cpu_count()])
        ncpus = 2

    # process the data either in parallel or serially
    if ncpus > 1:
        # prepare arguments for starmap
        items = []
        for i in range(len(con_files)):
            items.append((con_files[i], mag_files[i], dop_files[i], aia_files[i], mu_thresh, n_rings))
        # run in parellel
        print(">>> Processing %s epochs with %s processes..." % (len(con_files), ncpus))
        t0 = time.time()
        pids = []
        with get_context("spawn").Pool(ncpus, maxtasksperchild=2) as pool:
            # get PIDs of workers
            for child in mp.active_children():
                pids.append(child.pid)

            # run the analysis
            results = pool.starmap(process_data_set_parallel, items, chunksize=8)

        # find the output data sets
        datadir = str(root / "data") + "/"
        outfiles1 = glob.glob(datadir + "rv_full_disk_*")
        outfiles2 = glob.glob(datadir + "rv_regions_*")

        # stitch them together on the main process
        stitch_output_files(datadir + "rv_full_disk.csv", outfiles1, delete=True)
        stitch_output_files(datadir + "rv_regions.csv", outfiles2, delete=True)

        print("Parallel: --- %s seconds ---" % (time.time() - t0))
    else:
        # run serially
        print(">>> Processing %s epochs on a single process" % len(con_files))
        t0 = time.time()
        for i in range(len(con_files)):
            process_data_set(con_files[i], mag_files[i], dop_files[i], aia_files[i],
                             mu_thresh=mu_thresh, n_rings=n_rings, plot=False)
        print("Serial: --- %s seconds ---" % (time.time() - t0))
    return None

if __name__ == "__main__":
    main()
