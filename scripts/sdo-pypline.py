import numpy as np
import matplotlib.pyplot as plt
import pdb
import glob
import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning)

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# bring functions into scope
from sdo_pypline.sdo_image import *
from sdo_pypline.sdo_vels import *

# function to glob the input data
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

    # set mu threshold
    mu_thresh = 0.1

    # loop over files
    for i in range(len(con_files)):
        # make SDOImage instances
        con = SDOImage(con_files[i])
        mag = SDOImage(mag_files[i])
        dop = SDOImage(dop_files[i])
        aia = SDOImage(aia_files[i])

        # interpolate aia image onto hmi image scale
        aia.rescale_to_hmi(con)

        # correct magnetogram for foreshortening
        mag.correct_magnetogram()

        # calculate differential rotation & observer velocity
        dop.calc_vrot_vobs()

        # calculate limb darkening/brightening in continuum map and filtergram
        con.calc_limb_darkening()
        aia.calc_limb_darkening()

        # set values to nan for mu less than mu_thresh
        con.mask_low_mu(mu_thresh)
        dop.mask_low_mu(mu_thresh)
        mag.mask_low_mu(mu_thresh)
        aia.mask_low_mu(mu_thresh)

        # identify regions for thresholding
        mask = SunMask(con, mag, dop, aia)

        # compute velocities
        vels = calc_velocities(con, mag, dop, aia, mask)

        # compute velocities by region
        vels1 = calc_velocities(con, mag, dop, aia, mask, region=1)
        vels2 = calc_velocities(con, mag, dop, aia, mask, region=2)
        vels3 = calc_velocities(con, mag, dop, aia, mask, region=3)


        pdb.set_trace()

        # # compute velocities
        # vels = create_sun_mask(con, mag, dop, aia, mu_thresh=mu_thresh)
        # v_hat, v_phot, v_quiet, v_conv = vels

        # plot the data
        # mag.plot_image()
        # dop.plot_image()
        # con.plot_image()
        # aia.plot_image()

if __name__ == "__main__":
    main()
