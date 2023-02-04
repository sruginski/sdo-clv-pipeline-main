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

pl_color = "tab:purple" # colors[4]
nw_color = "tab:pink" # colors[6]
qs_color = "orange" # colors[1]
rp_color = "tab:red" # colors[3]
bp_color = "tab:blue" # colors[9]
pu_color = "sienna" # colors[5]
um_color = "tab:gray" # colors[7]

# get the sdo data to plot
con, mag, dop, aia = download_plot_data()

# reduce the data
print(">>> Processing and plotting data...")
con, mag, dop, aia, mask = reduce_sdo_images(con, mag, dop, aia)

# calculate velocities
v_hat = np.nansum(dop.v_corr * con.image * (con.mu >= 0.1)) / np.nansum(con.image * (con.mu >= 0.1))
v_quiet = np.nansum(dop.v_corr * con.image * mask.is_quiet_sun()) / np.nansum(con.image * mask.is_quiet_sun())
v_umb = np.nansum(dop.v_corr * con.image * mask.is_umbra()) / np.nansum(con.image * mask.is_umbra())
v_pen = np.nansum(dop.v_corr * con.image * mask.is_penumbra()) / np.nansum(con.image * mask.is_penumbra())
v_network = np.nansum(dop.v_corr * con.image * mask.is_network()) / np.nansum(con.image * mask.is_network())
v_plage = np.nansum(dop.v_corr * con.image * mask.is_plage()) / np.nansum(con.image * mask.is_plage())

v_conv = v_hat - v_quiet

# calculate fraction of light in each group
all_light = np.nansum(con.image * (con.mu >= 0.1))
umb_light = np.nansum(mask.is_umbra() * con.image)/all_light
pen_light = np.nansum(mask.is_penumbra() * con.image)/all_light
blu_pen_light = np.nansum(mask.is_blue_penumbra() * con.image)/all_light
red_pen_light = np.nansum(mask.is_red_penumbra() * con.image)/all_light
quiet_light = np.nansum(mask.is_quiet_sun() * con.image)/all_light
network_light = np.nansum(mask.is_network() * con.image)/all_light
plage_light = np.nansum(mask.is_plage() * con.image)/all_light

v_quiet2 = (v_quiet * umb_light + v_quiet * pen_light + v_quiet * network_light + v_quiet * plage_light)/(umb_light + pen_light + network_light + plage_light)

v_conv_umb = umb_light * (v_umb - v_quiet)
v_conv_pen = pen_light * (v_pen - v_quiet)
v_conv_network = network_light * (v_network - v_quiet)
v_conv_plage = plage_light * (v_plage - v_quiet)

v_conv_umb + v_conv_pen + v_conv_network + v_conv_plage

pdb.set_trace()
