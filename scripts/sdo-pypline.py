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

# actually do things
def main():
    # find the data
    con_files = glob.glob("/Users/michael/Desktop/sdo_data/*hmi.ic*.fits")
    mag_files = glob.glob("/Users/michael/Desktop/sdo_data/*hmi.m*.fits")
    dop_files = glob.glob("/Users/michael/Desktop/sdo_data/*hmi.v*.fits")
    aia_files = glob.glob("/Users/michael/Desktop/sdo_data/*aia*.fits")

    # make HDO_Image instances
    con_img1 = HMI_Image(con_files[0])
    mag_img1 = HMI_Image(mag_files[0])
    dop_img1 = HMI_Image(dop_files[0])
    aia_img1 = AIA_Image(aia_files[0])

    # mask low mus
    # con_img1.mask_low_mu(0.15)
    # mag_img1.mask_low_mu(0.15)
    # dop_img1.mask_low_mu(0.15)

    # correct magnetogram for foreshortening
    mag_img1.correct_magnetogram()

    # correct dopplergram for differential rotation & observer velocity
    dop_img1.correct_dopplergram()

    # correct limb darkening in continuum map and filtergram
    con_img1.correct_limb_darkening()
    # aia_img1.correct_limb_darkening()

    # interpolate aia image onto hmi image scale
    aia_img1.rescale_to_hmi(con_img1)

    # plot the data
    mag_img1.plot_image()
    dop_img1.plot_image()
    con_img1.plot_image()
    aia_img1.plot_image()

if __name__ == "__main__":
    main()
