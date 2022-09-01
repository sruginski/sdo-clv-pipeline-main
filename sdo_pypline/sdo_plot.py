import numpy as np
from .sdo_image import *

import pdb, warnings
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy import ndimage
from scipy.optimize import curve_fit
from reproject import reproject_interp
from astropy.wcs import FITSFixedWarning
from astropy.io.fits.verify import VerifyWarning


def plot_image(sdo_image, outdir=None):
    assert outdir is not None
    if sdo_image.is_magnetogram():
        # get cmap
        cmap = plt.get_cmap("RdYlBu").copy()
        cmap.set_bad(color="black")

        # plot the sun
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.imshow(sdo_image.image, cmap=cmap, origin="lower", vmin=-4200, vmax=4200, interpolation=None)
        cb = fig.colorbar(im)
        cb.set_label(r"${\rm Magnetic\ Field\ Strength\ (G)}$")
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(r"${\rm Corrected\ HMI\ LOS\ Magnetic\ Field}$")
        ax1.text(2650, 50, sdo_image.date_obs, fontsize=8, c="white")
        ax1.grid(False)
        fig.savefig(outdir + "mag_" + sdo_image.date_obs + ".pdf", bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    elif sdo_image.is_dopplergram():
        # get cmap
        cmap = plt.get_cmap("seismic").copy()
        cmap.set_bad(color="black")

        # plot the sun
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.imshow(sdo_image.image - sdo_image.v_rot - sdo_image.v_obs, origin="lower", cmap=cmap, vmin=-2000, vmax=2000, interpolation=None)
        cb = fig.colorbar(im)
        cb.set_label(r"${\rm LOS\ Velocity\ (km/s)}$")
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(r"${\rm Corrected\ HMI\ LOS\ Dopplergram}$")
        ax1.text(2650, 50, sdo_image.date_obs, fontsize=8, c="white")
        ax1.grid(False)
        fig.savefig(outdir + "dop_" + sdo_image.date_obs + ".pdf", bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    elif sdo_image.is_continuum():
        # get cmap
        cmap = plt.get_cmap("afmhot").copy()
        cmap.set_bad(color="black")

        # plot the sun
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.imshow(sdo_image.iflat, cmap=cmap, origin="lower", interpolation=None)#, vmin=20000)
        cb = fig.colorbar(im)
        cb.set_label(r"${\rm Relative\ Intensity}$")
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(r"${\rm Flattened\ HMI\ Continuum}$")
        ax1.text(2650, 50, sdo_image.date_obs, fontsize=8, c="white")
        ax1.grid(False)
        fig.savefig(outdir + "con_" + sdo_image.date_obs + ".pdf", bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    elif sdo_image.is_filtergram():
        # get cmap
        cmap = plt.get_cmap("Purples_r").copy()
        cmap.set_bad(color="black")

        # plot the sun
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.imshow(sdo_image.iflat, cmap=cmap, origin="lower", interpolation=None)#, vmin=20000)
        cb = fig.colorbar(im)
        cb.set_label(r"${\rm Relative\ Intensity}$")
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(r"${\rm Flattened\ AIA\ 1700\AA\ Filtergram}$")
        ax1.text(2650, 50, sdo_image.date_obs, fontsize=8, c="white")
        ax1.grid(False)
        fig.savefig(outdir + "aia_" + sdo_image.date_obs + ".pdf", bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    else:
        return None


def plot_mask(mask, outdir=None):
    assert outdir is not None

    # get cmap
    cmap = colors.ListedColormap(["black", "saddlebrown", "orange", "yellow"])
    cmap.set_bad(color="black")
    norm = colors.BoundaryNorm([0, 1, 2, 3, 4], ncolors=cmap.N, clip=True)

    # plot the sun
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    im = ax1.imshow(mask.regions - 0.5, cmap=cmap, norm=norm, origin="lower", interpolation=None)
    cb = fig.colorbar(im, ticks=[0.5, 1.5, 2.5, 3.5])
    cb.ax.set_yticklabels(["Umbrae", "Penumbrae", "Quiet Sun", "Plage"])
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    ax1.set_title(r"${\rm Identified\ Regions}$")
    ax1.text(2650, 50, mask.date_obs, fontsize=8, c="white")
    ax1.grid(False)
    fig.savefig(outdir + "mask_" + mask.date_obs + ".pdf", bbox_inches="tight", dpi=500)
    plt.clf(); plt.close()

    return None
