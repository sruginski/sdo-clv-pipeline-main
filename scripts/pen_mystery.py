import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from astropy.nddata import Cutout2D
from astropy import units as u

from sdo_pypline.paths import root
from sdo_pypline.sdo_plot import *
from sdo_pypline.sdo_process import *
from sdo_pypline.sdo_download import *
from download_plot_data import download_plot_data

datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

def main():
    # get the data
    start = "2014/01/05"
    end = "2014/01/09"
    sample = 24
    con, mag, dop, aia = download_data(outdir=datadir + "fits/", start=start, end=end, sample=sample)


    pdb.set_trace()

    # find contiguous penumbra regions
    structure = ndimage.generate_binary_structure(2,2)
    labels, nlabels = ndimage.label(indp, structure=structure)


    plt.imshow(con.iflat/con.ld_coeffs[0] * (labels == 919.0) * (mag.image > 1000.0), cmap="afmhot", alpha=0.5, vmin=0.0, vmax=1.1)
    plt.imshow(con.iflat/con.ld_coeffs[0] * (labels == 919.0) * (mag.image < 1000.0), cmap="Purples", alpha=0.5, vmin=0.0, vmax=1.1)
    plt.imshow(dop.v_corr * (labels == 919.0), cmap="seismic", vmin=-2500, vmax=2500)


    plt.hist((con.iflat/con.ld_coeffs[0])[labels == 919.0], density=True, bins="auto", histtype="step", label="center")
    plt.hist((con.iflat/con.ld_coeffs[0])[labels == 679.0], density=True, bins="auto", histtype="step", label="limb")

    plt.hist(mag.image[labels == 919.0], density=True, bins="auto", histtype="step", label="center")
    plt.hist(mag.image[labels == 679.0], density=True, bins="auto", histtype="step", label="limb")

    plt.hist(dop.v_corr[labels == 919.0] * con.image[labels == 919.0] / np.nansum(con.image[labels == 919.0]), density=True, bins="auto", histtype="step", label="center")
    plt.hist(dop.v_corr[labels == 679.0], density=True, bins="auto", histtype="step", label="limb")



    plt.legend()
    plt.show()

    pdb.set_trace()


    return None


if __name__ == "__main__":
    main()

