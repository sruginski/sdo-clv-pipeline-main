import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

from sdo_pypline.paths import root

# use style
plt.style.use("my.mplstyle"); plt.ioff()

def calc_whole_disk_vel(df):

    return None

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# read in the data and sort by mjd
df_all = pd.read_csv(datadir + "velocities.csv")
df_all.sort_values(by=["mjd"], inplace=True)

df_pixels = pd.read_csv(datadir + "pixel_stats.csv")
df_pixels.sort_values(by=["mjd"], inplace=True)

df_pixels = pd.read_csv(datadir + "light_stats.csv")
df_pixels.sort_values(by=["mjd"], inplace=True)

pdb.set_trace()

whole_disk_pixel = results_pixel[0]
whole_disk_light = results_light[0]

mus_light = results_light[1:]
mus_vels = results_vel[1:]

vel_whole_disk = results_vel[0][7]

vel = 0
derp = len(np.unique(mask.regions[~np.isnan(mask.regions)]))
for i in range(len(mus_light)):
    mus_light_i = mus_light[i]
    for k in np.unique(mask.regions[~np.isnan(mask.regions)]):
        vel += (mus_vels[i*derp + int(k) - 1][7] * mus_light_i[int(k)+2])
