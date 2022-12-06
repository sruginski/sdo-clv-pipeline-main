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
from download_plot_data import download_plot_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def plot_ld(img):
    # params for LD
    mu_lim=0.1
    num_mu=25

    # get average intensity in evenly spaced rings
    mu_edge = np.linspace(1.0, mu_lim, num=num_mu)
    mu_avgs = (mu_edge[1:] + mu_edge[0:-1]) / 2.0
    mu_step = np.abs(np.diff(mu_edge)[0])
    avg_int = np.zeros(len(mu_edge)-1)
    for i in range(len(avg_int)):
        # find indices in ring that aren't nan
        inds = (img.mu > mu_edge[i+1]) & (img.mu <= mu_edge[i]) & (~np.isnan(img.image))

        # mask section that are big outliers
        ints = img.image[inds]
        plt.scatter(mu_avgs[i], np.nanmean(ints), c="k")

        ints[np.abs(ints - np.mean(ints)) >= (2.0 * np.std(ints))] = np.nan
        avg_int[i] = np.nanmean(ints)

    plt.scatter(mu_avgs, avg_int, c="tab:blue")

    # set the initial guess parameters for optimization
    if img.is_continuum():
        p0 = [59000.0, 0.38, 0.23]
    else:
        p0 = [1000, 0.9, -0.25]

    # do the fit and divide out the LD profile
    popt, pcov = curve_fit(quad_darkening, mu_avgs, avg_int, p0=p0)
    ld_coeffs = popt
    ld_fit = popt[0] * quad_darkening_two(mu_avgs, *popt[1:])

    plt.plot(mu_avgs, ld_fit)
    plt.show()
    return None

def main():
    # get the sdo data to plot
    con, mag, dop, aia = download_plot_data()

    # get SDO image instances
    con = SDOImage(con)
    aia = SDOImage(aia)
    con.calc_geometry()
    aia.rescale_to_hmi(con)

    pdb.set_trace()

    plot_ld(con)
    plot_ld(aia)



if __name__ == "__main__":
    main()


