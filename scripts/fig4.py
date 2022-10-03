import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from scipy import ndimage
from skimage.measure import regionprops
from sdo_pypline.paths import root
from sdo_pypline.sdo_process import *
from sdo_pypline.sdo_download import *
from download_plot_data import download_plot_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def get_areas_for_set(con_file, mag_file, dop_file, aia_file):
    # reduce the data
    print(">>> Processing dataset...")
    con, mag, dop, aia, mask = reduce_sdo_images(con_file, mag_file, dop_file, aia_file, mu_thresh=0.1)

    # step through the region id process for AIA
    regions = np.zeros(np.shape(con.image))

    # calculate intensity thresholds for HMI
    con_thresh1 = 0.89 * np.nansum(con.iflat * mask.w_quiet)/np.nansum(mask.w_quiet)
    con_thresh2 = 0.45 * np.nansum(con.iflat * mask.w_quiet)/np.nansum(mask.w_quiet)

    # get indices for penumbrae, umbrae, quiet sun
    ind1 = con.iflat <= con_thresh2
    ind2 = (con.iflat <= con_thresh1) & (con.iflat > con_thresh2)
    ind3 = (con.iflat > con_thresh1) & mask.w_quiet

    # calculate intensity thresholds for AIA
    weights = mask.w_active * (~ind1) * (~ind2)
    aia_thresh = np.nansum(aia.iflat * weights)/np.nansum(weights)

    # get indices for bright regions (plage/faculae + network)
    ind4a = (con.iflat > con_thresh1) & mask.w_active
    ind4b = (aia.iflat > aia_thresh) & (~ind1) & (~ind2)
    ind4 = ind4a | ind4b

    # set mask indices
    regions[ind1] = 1 # umbrae
    regions[ind2] = 2 # penumbrae
    regions[ind3] = 3 # quiet sun
    regions[ind4] = 4 # bright areas (will separate into plage/faculae+network)

    # get number of pixels for area conversion to microhemispheres
    npix = np.sum(con.mu > 0.0)

    # label unique contiguous bright regions and calculate their sizes
    binary_img = regions == 4
    structure = ndimage.generate_binary_structure(2,2)
    labels, nlabels = ndimage.label(binary_img, structure=structure)

    # get labeled region properties
    rprops = regionprops(labels)
    areas = np.array([rprop.area for rprop in rprops]) / npix
    perimeters = np.array([rprop.perimeter for rprop in rprops])

    # get the latitude of the centroid
    centroids = np.zeros(len(areas))
    for i in range(len(centroids)):
        cent_row, cent_col = rprops[i].centroid
        cent_row = int(cent_row)
        cent_col = int(cent_col)
        centroids[i] = dop.phi[cent_row, cent_col]
    return areas, perimeters, centroids

def main():
    # see if the areas have already been saved
    areas_file = datadir + "areas.txt"
    perimeters_file = datadir + "perimeters.txt"
    centroids_file = datadir + "centroids.txt"
    if (not exists(areas_file)) | (not exists(perimeters_file)):
        # get the data to process (1 week at one day cadence)
        con_files = glob.glob(datadir + "fits/*cont*.fits")
        if len(con_files) < 30:
            start = "2014/01/01"
            end = "2014/01/31"
            sample = 24
            files = download_data(outdir=datadir+"fits/", start=start, end=end, sample=24, overwrite=False)
            con_files, mag_files, dop_files, aia_files = find_data(datadir+"fits/")
        else:
            con_files, mag_files, dop_files, aia_files = find_data(datadir+"fits/")

        # get the areas
        areas = []
        perimeters = []
        centroids = []
        for i in range(len(con_files)):
            if i > 10:
                break
            area, perimeter, centroid = get_areas_for_set(con_files[i], mag_files[i], dop_files[i], aia_files[i])
            areas.append(area)
            perimeters.append(perimeter)
            centroids.append(centroid)

        # write the areas to disk
        areas = np.concatenate(areas)
        perimeters = np.concatenate(perimeters)
        centroids = np.concatenate(centroids)
        np.savetxt(areas_file, areas)
        np.savetxt(perimeters_file, perimeters)
        np.savetxt(centroids_file, centroids)
    else:
        areas = np.loadtxt(areas_file)
        perimeters = np.loadtxt(perimeters_file)
        centroids = np.loadtxt(centroids_file)

    # calculate ratios
    ratios = perimeters/areas
    mask = (perimeters >= 2)

    pdb.set_trace()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(areas * 1e6, np.abs(np.sin(centroids - np.pi/2)))
    ax1.set_xscale("log")
    plt.show()

    # get the bins to use
    # bins = np.logspace(np.log10(np.min(perimeters[mask])), np.log10(np.max(perimeters)), 50)
    # plt.hist(perimeters[mask], bins=bins, histtype="step")
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xlabel("region perimeter (pixels)")
    # plt.ylabel("count")
    # plt.savefig("/Users/michael/Desktop/perimeter_dist.pdf")
    # plt.clf(); plt.close()

    # bins = np.logspace(np.log10(np.min(areas[mask])), np.log10(np.max(areas)), 50)
    # plt.hist(areas[mask], bins=bins, histtype="step")
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xlabel("region area (pixels)")
    # plt.ylabel("count")
    # plt.savefig("/Users/michael/Desktop/area_dist.pdf")
    # plt.clf(); plt.close()

    # bins = np.linspace(np.min(ratios[mask]), np.max(ratios), 50)
    # plt.hist(ratios[mask & (areas > 25)], bins=bins, histtype="step")
    # plt.yscale("log")
    # plt.xlabel("perimeter/area")
    # plt.ylabel("count")
    # plt.savefig("/Users/michael/Desktop/ratio_dist.pdf")
    # plt.clf(); plt.close()

    idx1 = np.random.choice(len(ratios[mask]), int(np.floor(len(ratios[mask])/1)))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axhline(30, ls="--", c="k", alpha=0.75)
    ax1.axvline(0.5, ls=":", c="k", alpha=0.75)
    ax1.scatter(ratios[mask][idx1], areas[mask][idx1], s=1, c="tab:blue", alpha=0.75, rasterized=True)
    # plt.scatter(ratios[mask][idx1], areas[mask][idx1], s=1, c=np.log10(perimeters[mask][idx1]), rasterized=True)
    # clb = plt.colorbar()
    # clb.set_label("log10(perimeter [pixels])")
    ax1.set_xlabel(r"${\rm perimeter/area}$")
    ax1.set_ylabel(r"${\rm region\ area\ (pixels)}$")
    ax1.set_yscale("log")
    # plt.legend()
    fig.savefig(plotdir + "fig4.pdf", dpi=200)
    plt.clf(); plt.close()

    return None

if __name__ == "__main__":
    main()
