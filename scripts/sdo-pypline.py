import numpy as np
import sunpy as sp
from sunpy.net import Fido, attrs as a
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb
import glob
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# bring functions into scope
from sdo_pypline.sdo_io import *
from sdo_pypline.sdo_image import *
# from .sdo_pypline.hmi_prep import *
# from .sdo_pypline.geometry import *

def find_sdo_data(indir):
    # find the data
    con_files = glob.glob(indir + "*hmi.ic*.fits")
    mag_files = glob.glob(indir + "*hmi.m*.fits")
    dop_files = glob.glob(indir + "*hmi.v*.fits")
    aia_files = glob.glob(indir + "*aia*.fits")
    return con_files, mag_files, dop_files, aia_files

# actually do things
def main():
    # find the data
    indir = "/Users/michael/Desktop/sdo_data/"
    con_files, mag_files, dop_files, aia_files = find_sdo_data(indir)

    # check the lengths
    assert (len(con_files) == len(mag_files) == len(dop_files) == len(aia_files))

    # TODO sort the files by date

    # loop over files
    for i in range(len(con_files)):
        # make SDOImage instances
        con_img1 = SDOImage(con_files[i])
        # mag_img1 = SDOImage(mag_files[i])
        # dop_img1 = SDOImage(dop_files[i])
        aia_img1 = SDOImage(aia_files[i])

        # mask low mus
        # con_img1.mask_low_mu(0.15)
        # mag_img1.mask_low_mu(0.15)
        # dop_img1.mask_low_mu(0.15)

        # correct magnetogram for foreshortening
        # mag_img1.correct_magnetogram()

        # correct dopplergram for differential rotation & observer velocity
        # dop_img1.correct_dopplergram()

        # correct limb darkening/brightening in continuum map and filtergram
        con_img1.correct_limb_darkening()
        aia_img1.correct_limb_darkening()

        # interpolate aia image onto hmi image scale
        aia_img1.rescale_to_hmi(con_img1)

        # plot the data
        # mag_img1.plot_image()
        # dop_img1.plot_image()
        # con_img1.plot_image()
        aia_img1.plot_image()

if __name__ == "__main__":
    main()
