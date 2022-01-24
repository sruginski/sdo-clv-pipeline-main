#==============================================================================
# Author: Michael Palumbo
# Date: May 2020
# Purpose:
#==============================================================================
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
    con_files = glob.glob("/Users/michael/Desktop/sdo_data/*con*.fits")
    mag_files = glob.glob("/Users/michael/Desktop/sdo_data/*mag*.fits")
    dop_files = glob.glob("/Users/michael/Desktop/sdo_data/*Dop*.fits")

    # make SDOImage instance
    con_img1 = HMI_Image(con_files[0])
    mag_img1 = HMI_Image(mag_files[0])
    dop_img1 = HMI_Image(dop_files[0])

    dop_img1.calc_mu_grid()

    # img = read_data(dat[0])
    # hdr = read_header(dat[0])



    pdb.set_trace()

    # # get grid of mus
    # int_dict = get_geom_params(hdr)
    # int_mu = get_mu_grid(int_dict)

    # # mask low mus
    # img[~(int_mu > 0.15)] = np.nan

    # # get cmap
    # cmap = plt.get_cmap("afmhot")
    # cmap.set_bad(color="white")

    # # plot the sun
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # im = ax1.imshow(img, cmap=cmap, origin="lower", vmin=20000)
    # # cb = fig.colorbar(im)
    # ax1.xaxis.set_visible(False)
    # ax1.yaxis.set_visible(False)
    # ax1.set_title(r"${\rm HMI\ Continuum\ Intensity}$")
    # ax1.text(2750, 50, hdr["DATE-OBS"], fontsize=8)
    # ax1.grid(False)
    # fig.savefig("/Users/michael/Desktop/test.pdf", bbox_inches="tight", dpi=500)
    # plt.clf(); plt.close()


    # # read in the data
    # dat = glob.glob("/Users/michael/Desktop/*mag*.fits")
    # img = read_data(dat[1])
    # hdr = read_header(dat[1])

    # # get grid of mus
    # int_dict = get_geom_params(hdr)
    # int_mu = get_mu_grid(int_dict)

    # # mask low mus
    # img[~(int_mu > 0.15)] = np.nan

    # # get cmap
    # cmap = plt.get_cmap("RdYlBu")
    # cmap.set_bad(color="white")

    # # plot the sun
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # im = ax1.imshow(img, cmap=cmap, origin="lower", vmin=-4200, vmax=4200)
    # # cb = fig.colorbar(im)
    # ax1.xaxis.set_visible(False)
    # ax1.yaxis.set_visible(False)
    # ax1.set_title(r"${\rm HMI\ LOS\ Magnetic\ Field\ Strength}$")
    # ax1.text(2750, 50, hdr["DATE-OBS"], fontsize=8)
    # ax1.grid(False)
    # fig.savefig("/Users/michael/Desktop/mag.pdf", bbox_inches="tight", dpi=500)
    # plt.clf(); plt.close()



if __name__ == "__main__":
    main()
