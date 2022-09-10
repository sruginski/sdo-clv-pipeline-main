import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import re, pdb, csv, glob, time, argparse
from astropy.time import Time
from os.path import exists, split, isdir, getsize

# bring functions into scope
from sdo_pypline.paths import root
from sdo_pypline.sdo_io import *
from sdo_pypline.sdo_vels import *
from sdo_pypline.sdo_image import *
from sdo_pypline.sdo_process import *

import multiprocessing as mp

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

# actually do things
def main():
    # initialize argparser
    parser = argparse.ArgumentParser(description="Analyze SDO data")
    parser.add_argument("--clobber", action="store_true", default=False)

    # parse the command line arguments
    args = parser.parse_args()
    clobber = args.clobber

    # define sdo_data directories
    indir = "/Users/michael/Desktop/sdo_data/"
    if not isdir(indir):
        indir = "/storage/home/mlp95/scratch/sdo_data/"

    # sort out input/output data files
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
        ncpus = mp.cpu_count()

    # process the data either in parallel or serially
    if ncpus > 1:
        # prepare arguments for starmap
        items = []
        for i in range(len(con_files)):
            items.append((con_files[i], mag_files[i], dop_files[i], aia_files[i], mu_thresh, n_rings))

        # figure out chunksize
        if len(items) <= ncpus:
            chunksize = 1
        else:
            chunksize = int(np.ceil(len(items)/ncpus))

        # run in parallel
        print(">>> About to parallel process with %s processes" % ncpus)
        t0 = time.time()
        with mp.Pool(ncpus) as pool:
            pool.starmap(process_data_set, items, chunksize=chunksize)
        print("Parallel: --- %s seconds ---" % (time.time() - t0))
    else:
        # run serially
        print(">>> About to serial process")
        t0 = time.time()
        for i in range(len(con_files)):
            process_data_set(con_files[i], mag_files[i], dop_files[i], aia_files[i],
                             mu_thresh=mu_thresh, n_rings=n_rings, plot=False)
        print("Serial: --- %s seconds ---" % (time.time() - t0))

if __name__ == "__main__":
    main()
