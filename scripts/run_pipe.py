import numpy as np
import matplotlib.pyplot as plt
import os, pdb, glob, time, argparse
from os.path import exists, split, isdir, getsize

# bring functions into scope
from sdo_pypline.paths import root
from sdo_pypline.sdo_io import *
from sdo_pypline.sdo_process import *

# multiprocessing imports
from multiprocessing import get_context
import multiprocessing as mp

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def get_parser_args():
    # initialize argparser
    parser = argparse.ArgumentParser(description="Analyze SDO data")
    parser.add_argument("--clobber", action="store_true", default=False)
    parser.add_argument("--globexp", type=str, default="")

    # parse the command line arguments
    args = parser.parse_args()
    clobber = args.clobber
    globexp = args.globexp
    return clobber, globexp

def main():
    # define sdo_data directories
    indir = "/Users/michael/Desktop/sdo-pypline/data/fits/"
    if not isdir(indir):
        indir = "/storage/home/mlp95/scratch/sdo_data/"

    # make raw data dir if it does not exist
    if not isdir(str(root / "data") + "/"):
        os.mkdir(str(root / "data") + "/")

    # sort out input/output data files
    clobber, globexp = get_parser_args()
    globdir = globexp.replace("*","")
    files = organize_IO(indir, clobber=clobber, globexp=globexp)
    con_files, mag_files, dop_files, aia_files = files

    # get output datadir
    datadir = str(root / "data") + "/" + globdir + "/"
    if not isdir(datadir):
        os.mkdir(datadir)

    # set mu threshold, number of mu rings
    n_rings = 10
    mu_thresh = 0.1
    plot = False

    # get number of cpus
    try:
        from os import sched_getaffinity
        print(">>> OS claims %s CPUs are available..." % len(sched_getaffinity(0)))
        ncpus = len(sched_getaffinity(0)) - 1
    except:
        # ncpus = np.min([len(con_files), mp.cpu_count()])
        ncpus = 1

    ncpus = 1

    # process the data either in parallel or serially
    if ncpus > 1:
        # make tmp directory
        tmpdir = datadir + "tmp/"
        if not isdir(tmpdir):
            os.mkdir(tmpdir)

        # prepare arguments for starmap
        items = []
        for i in range(len(con_files)):
            items.append((con_files[i], mag_files[i], dop_files[i], aia_files[i], mu_thresh, n_rings, datadir))

        # run in parellel
        print(">>> Processing %s epochs with %s processes..." % (len(con_files), ncpus))
        t0 = time.time()
        pids = []
        with get_context("spawn").Pool(ncpus, maxtasksperchild=4) as pool:
            # get PIDs of workers
            for child in mp.active_children():
                pids.append(child.pid)

            # run the analysis
            pool.starmap(process_data_set_parallel, items, chunksize=4)

        # find the output data sets
        outfiles1 = glob.glob(tmpdir + "thresholds_*")
        outfiles2 = glob.glob(tmpdir + "pixel_stats_*")
        outfiles3 = glob.glob(tmpdir + "light_stats_*")
        outfiles4 = glob.glob(tmpdir + "velocities_*")
        outfiles5 = glob.glob(tmpdir + "mag_stats_*")
        outfiles6 = glob.glob(tmpdir + "unweighted_velocities_*")
        outfiles7 = glob.glob(tmpdir + "intensities_*")
        outfiles8 = glob.glob(tmpdir + "average_velocities_*")

        # stitch them together on the main process
        delete = False
        stitch_output_files(datadir + "thresholds.csv", outfiles1, delete=delete)
        stitch_output_files(datadir + "pixel_stats.csv", outfiles2, delete=delete)
        stitch_output_files(datadir + "light_stats.csv", outfiles3, delete=delete)
        stitch_output_files(datadir + "velocities.csv", outfiles4, delete=delete)
        stitch_output_files(datadir + "mag_stats.csv", outfiles5, delete=delete)
        stitch_output_files(datadir + "unweighted_velocities.csv", outfiles6, delete=delete)
        stitch_output_files(datadir + "intensities.csv", outfiles7, delete=delete)
        stitch_output_files(datadir + "average_velocities.csv", outfiles8, delete=delete)

        # print run time
        print("Parallel: --- %s seconds ---" % (time.time() - t0))
    else:
        # run serially
        print(">>> Processing %s epochs on a single process" % len(con_files))
        t0 = time.time()
        for i in range(len(con_files)):
            process_data_set(con_files[i], mag_files[i], dop_files[i], aia_files[i],
                             mu_thresh=mu_thresh, n_rings=n_rings, datadir=datadir)

        # print run time
        print("Serial: --- %s seconds ---" % (time.time() - t0))
    return None

if __name__ == "__main__":
    main()
