import numpy as np
import pandas as pd
import os, pdb, glob, time, argparse
from os.path import exists, split, isdir, getsize

# bring functions into scope
from sdo_pypline.paths import root
from sdo_pypline.sdo_io import *
from sdo_pypline.sdo_process import *

def main():
    # figure out data directories
    files = []
    datadir = str(root / "data") + "/"

    # find the files to combine
    for file in os.listdir(datadir):
        file = datadir + file
        if not isdir(file):
            continue
        else:
            for f in glob.glob(file + "/" + "*.csv"):
                files.append(f)

    # names for output files
    fname1 = datadir + "intensities.csv"
    fname2 = datadir + "pixel_stats.csv"
    fname3 = datadir + "light_stats.csv"
    fname4 = datadir + "velocities.csv"
    fname5 = datadir + "mag_stats.csv"
    fname6 = datadir + "velocities_unweighted.csv"

    # headers for output files
    header1 = ["mjd", "aia_thresh", "a_aia", "b_aia", "c_aia", "hmi_thresh1", "hmi_thresh2", "a_hmi", "b_hmi", "c_hmi"]
    header2 = ["mjd", "lo_mu", "hi_mu", "ffactor", "umb_frac", "blu_pen_frac", "red_pen_frac", "quiet_frac", "network_frac", "plage_frac"]
    header3 = ["mjd", "lo_mu", "hi_mu", "umb_frac", "blu_pen_frac", "red_pen_frac", "quiet_frac", "network_frac", "plage_frac"]
    header4 = ["mjd", "region", "lo_mu", "hi_mu", "v_hat", "v_phot", "v_quiet", "v_conv"]
    header5 = ["mjd", "region", "lo_mu", "hi_mu", "mag_avg", "mag_std", "mag_unsigned"]
    header6 = ["mjd", "region", "lo_mu", "hi_mu", "v_hat", "v_phot", "v_quiet", "v_conv"]

    # delete old files if they exists
    fileset = (fname1, fname2, fname3, fname4, fname5, fname6)
    clean_output_directory(*fileset)

    # create the files with headers
    create_file(fname1, header1)
    create_file(fname2, header2)
    create_file(fname3, header3)
    create_file(fname4, header4)
    create_file(fname5, header5)
    create_file(fname6, header6)

    # now loop over files to combine
    for file in fileset:
        # find files to merge into output
        in_fname = []
        for f in files:
            if os.path.split(f)[-1] == os.path.split(file)[-1]:
                in_fname.append(f)

        # read in each file sequentially and append its data to output file
        for f in in_fname:
            data = pd.read_csv(f)
            data.to_csv(file, mode="a", index=False, header=False)

    return None


if __name__ == "__main__":
    main()

