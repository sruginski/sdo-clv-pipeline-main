import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from scipy import ndimage
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

    # label unique contiguous bright regions and calculate their sizes
    binary_img = regions == 4
    structure = ndimage.generate_binary_structure(2,2)
    labels, nlabels = ndimage.label(binary_img, structure=structure)
    areas = ndimage.sum(binary_img, labels, range(nlabels+1))
    return areas

def main():
    # see if the areas have already been saved
    areas_file = datadir + "areas.txt"
    if not exists(areas_file):
        # get the data to process (1 week at one day cadence)
        start = "2014/01/01"
        end = "2014/01/21"
        sample = 24
        files = download_data(outdir=datadir+"fits/", start=start, end=end, sample=24, overwrite=False)
        con_files, mag_files, dop_files, aia_files = find_data(datadir+"fits/")

        # get the areas
        areas = []
        for i in range(len(con_files)):
            print(i)
            areas.append(get_areas_for_set(con_files[i], mag_files[i], dop_files[i], aia_files[i]))

        # write the areas to disk
        areas = np.concatenate(areas)
        np.savetxt(areas_file, areas)
    else:
        areas = np.loadtxt(areas_file)

    pdb.set_trace()

    # only use areas > 1 pix
    areas = areas[areas > 1]

    # get the bins to use
    linbins = np.log10(np.linspace(4, 10, 5))
    linbins = []
    logbins = np.logspace(np.log10(10),np.log10(np.max(areas)),100)
    bins = np.concatenate((linbins, logbins))

    # get the distribution
    ahist, bin_edges = np.histogram(areas, bins=bins)
    bin_centers = (bin_edges[1:] + bin_edges[0:-1]) / 2
    indices = np.arange(1, len(ahist) + 1)

    # mask the nans and find where distribution goes to 0
    # TODO "polyfit may be poorly conditioned"
    peak_idx = np.argmax(ahist)
    zero_idx = np.argmax(bin_centers >= 1e4)

    # get data to fit
    mask = ahist[peak_idx:zero_idx] == 0.0
    xs = np.log10(bin_centers[peak_idx:zero_idx][~mask])
    ys = np.log10(ahist[peak_idx:zero_idx][~mask])
    ws = 1.0/np.sqrt(ahist[peak_idx:zero_idx][~mask])
    order = 1

    # fit the distribution
    # TODO should we weight?
    pfit1 = np.polyfit(xs, ys, order, w=ws)

    # evaluate the model and find where the fit cuts off the large-area tail
    model_ys = np.polyval(pfit1, np.log10(bin_centers))
    thresh_idx = np.argmax(model_ys < 0.0)
    area_thresh = bin_centers[thresh_idx]

    # make the plot
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(bin_centers[peak_idx])
    ax1.axvline(bin_centers[zero_idx])
    ax1.hist(areas, bins=bins, histtype="step")
    ax1.plot(bin_centers, 10**model_ys)
    ax1.set_xlim(1, ax1.get_xlim()[1])
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Region area (pixels)")
    ax1.set_ylabel("Count")

    pdb.set_trace()
    print("derp")
    print("derp")
    print("derp")
    print("derp")
    print("derp")

    return None

if __name__ == "__main__":
    main()
