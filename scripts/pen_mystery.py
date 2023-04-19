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
    # # get the data
    # start = "2014/01/05"
    # end = "2014/01/11"
    # sample = 24
    # con, mag, dop, aia = download_data(outdir=datadir + "fits/", start=start, end=end, sample=sample)

    # get file names
    con_file1 = datadir + "fits/hmi_ic_45s_2014_01_05_00_01_30_tai_continuum.fits"
    dop_file1 = datadir + "fits/hmi_v_45s_2014_01_05_00_01_30_tai_dopplergram.fits"
    mag_file1 = datadir + "fits/hmi_m_45s_2014_01_05_00_01_30_tai_magnetogram.fits"
    aia_file1 = datadir + "fits/aia_lev1_1700a_2014_01_05t00_00_30_74z_image_lev1.fits"

    con_file2 = datadir + "fits/hmi_ic_45s_2014_01_08_00_47_15_tai_continuum.fits"
    dop_file2 = datadir + "fits/hmi_v_45s_2014_01_08_00_47_15_tai_dopplergram.fits"
    mag_file2 = datadir + "fits/hmi_m_45s_2014_01_08_00_47_15_tai_magnetogram.fits"
    aia_file2 = datadir + "fits/aia_lev1_1700a_2014_01_08t00_00_30_71z_image_lev1.fits"

    con_file3 = datadir + "fits/hmi_ic_45s_2014_01_11_00_01_30_tai_continuum.fits"
    dop_file3 = datadir + "fits/hmi_v_45s_2014_01_11_00_01_30_tai_dopplergram.fits"
    mag_file3 = datadir + "fits/hmi_m_45s_2014_01_11_00_01_30_tai_magnetogram.fits"
    aia_file3 = datadir + "fits/aia_lev1_1700a_2014_01_11t00_00_30_71z_image_lev1.fits"

    con1, mag1, dop1, aia1, mask1 = reduce_sdo_images(con_file1, mag_file1, dop_file1, aia_file1)
    con2, mag2, dop2, aia2, mask2 = reduce_sdo_images(con_file2, mag_file2, dop_file2, aia_file2)
    con3, mag3, dop3, aia3, mask3 = reduce_sdo_images(con_file3, mag_file3, dop_file3, aia_file3)

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

