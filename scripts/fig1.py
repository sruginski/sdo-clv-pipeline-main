import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from sdo_pypline.paths import root
from sdo_pypline.sdo_plot import *
from sdo_pypline.sdo_process import *
from sdo_pypline.sdo_download import download_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def main():
    # see if data are already downloaded
    files = glob.glob(datadir + "*.fits")
    con = datadir + "hmi_ic_45s_2014_01_07_00_01_30_tai_continuum.fits"
    mag = datadir + "hmi_m_45s_2014_01_07_00_01_30_tai_magnetogram.fits"
    dop = datadir + "hmi_v_45s_2014_01_07_00_01_30_tai_dopplergram.fits"
    aia = datadir + "aia_lev1_1700a_2014_01_07t00_00_30_71z_image_lev1.fits"

    if any(map(lambda x: x not in files, (con, mag, dop, aia))):
        # download the data to plot
        print("\t >>> Downloading data")
        start = "2014/01/07"
        end = "2014/01/07"
        sample = 24
        con, mag, dop, aia = download_data(outdir=datadir + "fits/",
                                           start=start, end=end, sample=sample)

    # preprocess the data and plot it
    print("\t >>> Processing and plotting data...")
    process_data_set(con, mag, dop, aia, plot=True, vels=False)

    # find and rename the files
    con_file = glob.glob(plotdir + "*con*")
    mag_file = glob.glob(plotdir + "*mag*")
    dop_file = glob.glob(plotdir + "*dop*")
    aia_file = glob.glob(plotdir + "*aia*")
    mask_file = glob.glob(plotdir + "*mask*")

    os.rename(con_file[0], plotdir+"fig1a.pdf")
    os.rename(mag_file[0], plotdir+"fig1b.pdf")
    os.rename(dop_file[0], plotdir+"fig1c.pdf")
    os.rename(aia_file[0], plotdir+"fig1d.pdf")
    os.rename(mask_file[0], plotdir+"fig1e.pdf")
    return None

if __name__ == "__main__":
    main()


