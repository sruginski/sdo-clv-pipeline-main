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
from download_plot_data import download_plot_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def main():
    # get the sdo data to plot
    con, mag, dop, aia = download_plot_data()

    # reduce the data
    print(">>> Processing and plotting data...")
    con, mag, dop, aia, mask = reduce_sdo_images(con, mag, dop, aia)

    # plot them
    print("\t>>> Plotting HMI continuum...")
    plot_image(con, outdir=plotdir, fname="fig1a.pdf")

    print("\t>>> Plotting magnetogram...")
    plot_image(mag, outdir=plotdir, fname="fig1b.pdf")

    print("\t>>> Plotting dopplergram...")
    plot_image(dop, outdir=plotdir, fname="fig1c.pdf")

    print("\t>>> Plotting AIA continuum...")
    plot_image(aia, outdir=plotdir, fname="fig1d.pdf")

    print("\t>>> Plotting mask...")
    plot_mask(mask, outdir=plotdir, fname="fig3.pdf")
    return None

if __name__ == "__main__":
    main()


