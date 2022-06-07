import numpy as np
import matplotlib.pyplot as plt
import pdb
import csv
import glob
from astropy.time import Time
from os.path import exists, split, isdir

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

def write_vels(fname, mjd, vels):
    # create the file if it doesn't exist
    if (~exists(fname) & isdir(split(fname)[0])):
        with open(fname, "w") as f:
            writer = csv.writer(f)
            writer.writerow(["mjd", "v_hat", "v_phot", "v_quiet", "v_conv"])

    # write the vels
    with open(fname, "a") as f:
        writer = csv.writer(f)
        writer.writerow(np.concatenate(([mjd], [v for v in vels])))

    return None

def write_vels_by_region(fname, mjd, region, lo_mu, hi_mu, vels):
    # create the file if it doesn't exist
    if (~exists(fname) & isdir(split(fname)[0])):
        with open(fname, "w") as f:
            writer = csv.writer(f)
            writer.writerow(["mjd", "region", "lo_mu", "hi_mu", "v_hat", "v_phot", "v_quiet", "v_conv"])

    # write the vels
    with open(fname, "a") as f:
        writer = csv.writer(f)
        writer.writerow(np.concatenate(([mjd], [region], [lo_mu], [hi_mu], [v for v in vels])))

    return None

# actually do things
def main():
    # find the data
    indir = "/Users/michael/Desktop/sdo_data/"
    outdir = "/Users/michael/Desktop/"
    if not isdir(indir):
        indir = "/storage/home/mlp95/scratch/sdo_data/"
        outdir = "/storage/home/mlp95/work/sdo_output/"

    con_files, mag_files, dop_files, aia_files = find_sdo_data(indir)

    # check the lengths
    assert (len(con_files) == len(mag_files) == len(dop_files) == len(aia_files))

    # TODO sort the files by date

    # set mu threshold and number of mu_rings
    n_rings = 12
    mu_thresh = 0.1

    # grid of mu values to iterate over
    # TODO, evenly spaced in mu or theta?
    mu_grid = np.linspace(mu_thresh, 1.0, n_rings)

    # loop over files
    for i in range(len(con_files)):
        # make SDOImage instances
        con = SDOImage(con_files[i])
        mag = SDOImage(mag_files[i])
        dop = SDOImage(dop_files[i])
        aia = SDOImage(aia_files[i])

        # get MJD for observations
        mjd = Time(con.date_obs).mjd

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

        # compute velocities and write to disk
        vels = calc_velocities(con, mag, dop, aia, mask)
        write_vels(outdir + "rv_full_disk.csv", mjd, vels)

        # plot the data
        # mag.plot_image()
        # dop.plot_image()
        # con.plot_image()
        # aia.plot_image()

        # loop over mu annuli
        for i in range(n_rings-1):
            # mu values for annuli
            lo_mu=mu_grid[i]
            hi_mu=mu_grid[i+1]

            # loop over unique region identifiers
            for j in np.unique(mask.regions[~np.isnan(mask.regions)]):
                # compute velocity components in each mu annulus
                vels_reg = calc_velocities(con, mag, dop, aia, mask, region=j, lo_mu=lo_mu, hi_mu=hi_mu)

                # write to disk
                write_vels_by_region(outdir + "rv_regions.csv", mjd, j, lo_mu, hi_mu, vels_reg)

if __name__ == "__main__":
    main()
