import numpy as np
from .sdo_io import *
from .sdo_image import *

import pdb, warnings
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from astropy.io import fits
from astropy.wcs import WCS

from astropy.wcs import FITSFixedWarning
from astropy.io.fits.verify import VerifyWarning


def plot_image(sdo_image, outdir=None, fname=None):
    assert outdir is not None

    # get the WCS
    wcs = sdo_image.wcs

    # initialize the figure
    fig = plt.figure(figsize=(6.4, 4.8))
    ax1 = fig.add_subplot(111, projection=wcs)

    if sdo_image.is_magnetogram():
        # get cmap
        cmap = plt.get_cmap("RdYlBu").copy()
        cmap.set_bad(color="black")

        # plot the sun
        img = ax1.imshow(sdo_image.image, cmap=cmap, origin="lower", vmin=-4200, vmax=4200, interpolation=None)
        clb = fig.colorbar(img)
        clb.set_label(r"${\rm Magnetic\ Field\ Strength\ (G)}$")
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
        ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
        ax1.text(1400, 4000, sdo_image.date_obs, fontsize=10, c="white")
        ax1.grid(False)

        # figure out the filename
        if fname is None:
            fname = "mag_" + sdo_image.date_obs + ".pdf"
        fig.savefig(outdir + fname, bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    elif sdo_image.is_dopplergram():
        # get cmap
        cmap = plt.get_cmap("seismic").copy()
        cmap.set_bad(color="black")

        # plot the sun
        img = ax1.imshow(sdo_image.image - sdo_image.v_rot - sdo_image.v_obs, origin="lower", cmap=cmap, vmin=-2000, vmax=2000, interpolation=None)
        clb = fig.colorbar(img)
        clb.set_label(r"${\rm LOS\ Velocity\ (km/s)}$")
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
        ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
        ax1.text(1400, 4000, sdo_image.date_obs, fontsize=10, c="white")
        ax1.grid(False)

        # figure out the filename
        if fname is None:
            fname = "dop_" + sdo_image.date_obs + ".pdf"
        fig.savefig(outdir + fname, bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    elif sdo_image.is_continuum():
        # get cmap
        cmap = plt.get_cmap("afmhot").copy()
        cmap.set_bad(color="black")

        # plot the sun
        img = ax1.imshow(sdo_image.iflat, cmap=cmap, origin="lower", interpolation=None)#, vmin=20000)
        clb = fig.colorbar(img)
        clb.set_label(r"${\rm Relative\ Intensity}$")
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
        ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
        ax1.text(1400, 4000, sdo_image.date_obs, fontsize=10, c="white")
        ax1.grid(False)

        # figure out the filename
        if fname is None:
            fname = "con_" + sdo_image.date_obs + ".pdf"
        fig.savefig(outdir + fname, bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    elif sdo_image.is_filtergram():
        # get cmap
        cmap = plt.get_cmap("Purples_r").copy()
        cmap.set_bad(color="black")

        # plot the sun
        img = ax1.imshow(sdo_image.iflat, cmap=cmap, origin="lower", interpolation=None)#, vmin=20000)
        clb = fig.colorbar(img)
        clb.set_label(r"${\rm Relative\ Intensity}$")
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
        ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
        ax1.text(1400, 4000, sdo_image.date_obs, fontsize=10, c="white")
        ax1.grid(False)

        # figure out the filename
        if fname is None:
            fname = "aia_" + sdo_image.date_obs + ".pdf"
        fig.savefig(outdir + fname, bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None

    else:
        return None

def plot_mask(mask, outdir=None, fname=None):
    assert outdir is not None

    # get cmap
    cmap = colors.ListedColormap(["black", "saddlebrown", "orange", "yellow", "white"])
    cmap.set_bad(color="black")
    norm = colors.BoundaryNorm([0, 1, 2, 3, 4, 5], ncolors=cmap.N, clip=True)

    # get the WCS
    wcs = mask.wcs

    # plot the sun
    fig = plt.figure(figsize=(6.4, 4.8))
    ax1 = fig.add_subplot(111, projection=wcs)
    img = ax1.imshow(mask.regions - 0.5, cmap=cmap, norm=norm, origin="lower", interpolation=None)
    clb = fig.colorbar(img, ticks=[0.5, 1.5, 2.5, 3.5, 4.5])
    clb.ax.set_yticklabels([r"${\rm Umbra}$", r"${\rm Penumbra}$", r"${\rm Quiet\ Sun}$", r"${\rm Network}$", r"${\rm Plage}$"])
    ax1.invert_xaxis()
    ax1.invert_yaxis()
    ax1.set_xlabel(r"${\rm Helioprojective\ Longitude}$")
    ax1.set_ylabel(r"${\rm Helioprojective\ Latitude}$")
    ax1.text(1400, 4000, mask.date_obs, fontsize=10, c="white")
    ax1.grid(False)

    # figure out the filename
    if fname is None:
        fname = "mask_" + mask.date_obs + ".pdf"
    fig.savefig(outdir + fname, bbox_inches="tight", dpi=500)
    plt.clf(); plt.close()
    return None
