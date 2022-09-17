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

    con_files=con_files[0:4]
    mag_files=mag_files[0:4]
    dop_files=dop_files[0:4]
    aia_files=aia_files[0:4]

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
        ncpus = 4

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
            pool.starmap(process_data_set_parallel, items, chunksize=1)

        # find the output data sets
        datadir = str(root / "data") + "/"
        tmpdir = datadir + "tmp/"
        outfiles1 = glob.glob(tmpdir + "rv_full_disk_*")
        outfiles2 = glob.glob(tmpdir + "rv_regions_*")
        outfiles3 = glob.glob(tmpdir + "rv_mu_*")
        outfiles4 = glob.glob(tmpdir + "*aia*")
        outfiles5 = glob.glob(tmpdir + "*hmi*")
        outfiles6 = glob.glob(tmpdir + "*thresh*")

        # stitch them together on the main process
        stitch_output_files(datadir + "rv_full_disk.csv", outfiles1, delete=True)
        stitch_output_files(datadir + "rv_regions.csv", outfiles2, delete=True)
        stitch_output_files(datadir + "rv_mu.csv", outfiles3, delete=True)
        stitch_output_files(datadir + "aia_ld_params.csv", outfiles4, delete=True)
        stitch_output_files(datadir + "hmi_ld_params.csv", outfiles5, delete=True)
        stitch_output_files(datadir + "con_thresh.csv", outfiles6, delete=True)

        print("Parallel: --- %s seconds ---" % (time.time() - t0))
    else:
        # run serially
        print(">>> Processing %s epochs on a single process" % len(con_files))
        t0 = time.time()
        for i in range(len(con_files)):
            process_data_set(con_files[i], mag_files[i], dop_files[i],
                             aia_files[i], mu_thresh=mu_thresh, n_rings=n_rings)
        print("Serial: --- %s seconds ---" % (time.time() - t0))
    return None

if __name__ == "__main__":
    main()
