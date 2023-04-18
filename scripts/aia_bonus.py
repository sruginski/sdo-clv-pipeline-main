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
from sdo_pypline.sdo_download import *
from download_plot_data import download_plot_data

datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"


def main():
    pdb.set_trace()

    ind5a = ind5a.astype(float)
    ind5b = ind5b.astype(float)

    derp1 = np.zeros_like(ind5a)
    derp1[(ind5a == 0.0) & (ind5b == 1.0)] = 1.0

    derp2 = np.zeros_like(ind5a)
    derp2[(ind5a == 1.0) & (ind5b == 0.0)] = 1.0


    cmapa = plt.get_cmap("Oranges").copy()
    cmapa.set_bad(color="none")


    cmapb = plt.get_cmap("Purples").copy()
    cmapb.set_bad(color="none")

    plt.imshow(ind5a, cmap=cmapa)
    plt.imshow(ind5b, cmap=cmapb)
    plt.show()

    return None


