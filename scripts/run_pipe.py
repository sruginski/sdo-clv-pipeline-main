import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import re, pdb, csv, glob, argparse
from astropy.time import Time
from os.path import exists, split, isdir

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# bring functions into scope
from sdo_pypline.sdo_image import *
from sdo_pypline.sdo_vels import *

# actually do things
def main():
    # sort out directories
    indir = "/Users/michael/Desktop/sdo_data/"
    outdir = "/Users/michael/Desktop/"
    if not isdir(indir):
        indir = "/storage/home/mlp95/scratch/sdo_data/"
        outdir = "/storage/home/mlp95/work/sdo_output/"

    # set file names and truncate them if they exist
    fname1 = truncate_file(outdir + "rv_full_disk.csv")
    fname2 = truncate_file(outdir + "rv_mu.csv")
    fname3 = truncate_file(outdir + "rv_regions.csv")

    # find the input data
    con_files, mag_files, dop_files, aia_files = find_data(indir)

    # check the lengths
    assert (len(con_files) == len(mag_files) == len(dop_files) == len(aia_files))

    # set mu threshold and number of mu_rings
    n_rings = 10
    mu_thresh = 0.1

    # grid of mu values to iterate over
    # TODO, evenly spaced in mu or theta?
    mu_grid = np.linspace(mu_thresh, 1.0, n_rings)

    # decide whether to plot
    plot = False

    # loop over files
    for i in range(len(con_files)):
        # make SDOImage instances
        try:
            con = SDOImage(con_files[i])
            mag = SDOImage(mag_files[i])
            dop = SDOImage(dop_files[i])
            aia = SDOImage(aia_files[i])
        except OSError:
            print("\t >>> Invalid file, skipping " + get_date(con_files[i]).isoformat())
            continue

        # get MJD for observations and report status
        iso = Time(con.date_obs).iso
        mjd = Time(con.date_obs).mjd

        # interpolate aia image onto hmi image scale
        aia.rescale_to_hmi(con)

        # correct magnetogram for foreshortening
        mag.correct_magnetogram()

        # calculate differential rotation & observer velocity
        dop.calc_vrot_vobs()

        # calculate limb darkening/brightening in continuum map and filtergram
        try:
            con.calc_limb_darkening()
            aia.calc_limb_darkening()
        except RuntimeError:
            print("\t >>> Limb darkening fit failed, skipping " + iso)

        # set values to nan for mu less than mu_thresh
        con.mask_low_mu(mu_thresh)
        dop.mask_low_mu(mu_thresh)
        mag.mask_low_mu(mu_thresh)
        aia.mask_low_mu(mu_thresh)

        # identify regions for thresholding
        mask = SunMask(con, mag, dop, aia)

        # plot the data
        if plot:
            mag.plot_image(outdir=outdir)
            dop.plot_image(outdir=outdir)
            con.plot_image(outdir=outdir)
            aia.plot_image(outdir=outdir)
            mask.plot_image(outdir=outdir)

        # compute velocities and write to disk
        vels = calc_velocities(con, mag, dop, aia, mask)
        write_vels(fname1, mjd, mask.ff, mask.pen_frac,
                   mask.umb_frac, mask.quiet_frac,
                   mask.plage_frac, vels)

        # loop over mu annuli
        for i in range(n_rings-1):
            # mu values for annuli
            lo_mu=mu_grid[i]
            hi_mu=mu_grid[i+1]

            # compute velocity in mu annulus
            vels_reg = calc_velocities(con, mag, dop, aia, mask, lo_mu=lo_mu, hi_mu=hi_mu)

            # write to disk
            write_vels_by_region(fname2, mjd, 0, lo_mu, hi_mu, vels_reg)

            # loop over unique region identifiers
            for j in np.unique(mask.regions[~np.isnan(mask.regions)]):
                # compute velocity components in each mu annulus by region
                vels_reg = calc_velocities(con, mag, dop, aia, mask, region=j, lo_mu=lo_mu, hi_mu=hi_mu)

                # write to disk
                write_vels_by_region(fname3, mjd, j, lo_mu, hi_mu, vels_reg)

        # report status
        print("\t >>> Epoch " + iso + " run successfully!")

if __name__ == "__main__":
    main()
