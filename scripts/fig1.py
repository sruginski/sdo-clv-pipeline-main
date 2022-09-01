import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

# sort out paths
from sdo_pypline.paths import user as Paths
plotdir = str(Paths().figures) + "/"

# use style
plt.style.use("my.mplstyle"); plt.ioff()
