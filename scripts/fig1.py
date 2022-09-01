import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

from sdo_pypline.paths import root
from sdo_pypline.sdo_io import *
from sdo_pypline.sdo_image import *
from sdo_pypline.sdo_download import download_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

# download the data to plot
start = "2014/01/07"
end = "2014/01/07"
sample = 24
download_data(outdir=datadir, start=start, end=end, sample=sample)
