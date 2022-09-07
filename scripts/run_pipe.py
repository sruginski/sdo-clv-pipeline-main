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
    indir = "/Users/michael/work/sdo-pypline/data/"
    # outdir = "/Users/michael/Desktop/"
    if not isdir(indir):
        indir = "/storage/home/mlp95/scratch/sdo_data/"
        # outdir = "/storage/home/mlp95/work/sdo_output/"

    # sort out input/output data files
    con_files, mag_files, dop_files, aia_files = organize_input_output(indir, clobber=clobber)

    # set mu threshold, number of mu rings
    n_rings = 10
    mu_thresh = 0.1
    plot = False

    # loop over files
    for i in range(len(con_files)):
        # report status
        print("\t >>> Running epoch " + get_date(con_files[i]).isoformat())

        # analyze set of files
        process_data_set(con_files[i], mag_files[i], dop_files[i], aia_files[i],
                         mu_thresh=mu_thresh, n_rings=n_rings, plot=True)


if __name__ == "__main__":
    main()
