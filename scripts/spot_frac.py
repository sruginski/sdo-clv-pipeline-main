import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

from sdo_pypline.paths import root

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"
procdir = datadir + "processed/"


# read in by region
df_vels_full = pd.read_csv(procdir + "full_disk.csv")
plage = pd.read_csv(procdir + "plage.csv")
network = pd.read_csv(procdir + "network.csv")
quiet_sun = pd.read_csv(procdir + "quiet_sun.csv")
penumbrae = pd.read_csv(procdir + "penumbrae.csv")
red_penumbrae = pd.read_csv(procdir + "red_penumbrae.csv")
blu_penumbrae = pd.read_csv(procdir + "blu_penumbrae.csv")
umbrae = pd.read_csv(procdir + "umbrae.csv")

unq_mjds = np.unique(df_vels_full.mjd)

umb_fracs = []
pen_fracs = []
pla_fracs = []
net_fracs = []
for i, mjd in enumerate(unq_mjds):
    umb_fracs.append(np.sum(umbrae.pixel_frac[umbrae.mjd == unq_mjds[i]]))
    pen_fracs.append(np.sum(penumbrae.pixel_frac[penumbrae.mjd == unq_mjds[i]]))
    pla_fracs.append(np.sum(plage.pixel_frac[plage.mjd == unq_mjds[i]]))
    net_fracs.append(np.sum(network.pixel_frac[network.mjd == unq_mjds[i]]))

umb_fracs = np.array(umb_fracs)
pen_fracs = np.array(pen_fracs)
pla_fracs = np.array(pla_fracs)
net_fracs = np.array(net_fracs)

plt.scatter(unq_mjds, umb_fracs, s=0.5, label="Umbrae")
plt.scatter(unq_mjds, pen_fracs, s=0.5, label="Penmbrae")
plt.scatter(unq_mjds, umb_fracs + pen_fracs, s=0.5, label="Entire spot")
plt.scatter(unq_mjds, pla_fracs, s=0.5, label="Plage")
plt.scatter(unq_mjds, net_fracs, s=0.5, label="Network")
plt.xlabel("MJD")
plt.ylabel("Pixel fraction")
plt.legend()
plt.show()


pdb.set_trace()
