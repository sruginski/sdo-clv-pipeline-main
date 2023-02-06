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

# calculate whole-disk velocities
v_hat_whole = np.nansum(dop.v_corr * con.image * (con.mu >= 0.1)) / np.nansum(con.image * (con.mu >= 0.1))
v_quiet_whole = np.nansum(dop.v_corr * con.image * mask.is_quiet_sun()) / np.nansum(con.image * mask.is_quiet_sun())
v_umb_whole = np.nansum(dop.v_corr * con.image * mask.is_umbra()) / np.nansum(con.image * mask.is_umbra())
v_pen_whole = np.nansum(dop.v_corr * con.image * mask.is_penumbra()) / np.nansum(con.image * mask.is_penumbra())
v_network_whole = np.nansum(dop.v_corr * con.image * mask.is_network()) / np.nansum(con.image * mask.is_network())
v_plage_whole = np.nansum(dop.v_corr * con.image * mask.is_plage()) / np.nansum(con.image * mask.is_plage())

# get v_conv for whole disk
v_conv_whole = v_hat_whole - v_quiet_whole
v_conv_umb_whole = v_umb_whole - v_quiet_whole
v_conv_pen_whole = v_pen_whole - v_quiet_whole
v_conv_network_whole = v_network_whole - v_quiet_whole
v_conv_plage_whole = v_plage_whole - v_quiet_whole

# calculate fraction of light in each group
all_light = np.nansum(con.image * (con.mu >= 0.1))
umb_light = np.nansum(mask.is_umbra() * con.image)/all_light
pen_light = np.nansum(mask.is_penumbra() * con.image)/all_light
quiet_light = np.nansum(mask.is_quiet_sun() * con.image)/all_light
network_light = np.nansum(mask.is_network() * con.image)/all_light
plage_light = np.nansum(mask.is_plage() * con.image)/all_light

# combine v_conv in regions
v_conv_whole1 = umb_light * v_conv_umb_whole + pen_light * v_conv_pen_whole + network_light * v_conv_network_whole + plage_light * v_conv_plage_whole
v_conv_whole - v_conv_whole1

# now do it in two mu bins
mu_reg1 = (con.mu >= 0.5)
mu_reg2 = (con.mu >= 0.1) & (con.mu < 0.5)

# get velocities in bins
v_hat_reg1 = np.nansum(dop.v_corr * con.image * (con.mu >= 0.1) * mu_reg1) / np.nansum(con.image * (con.mu >= 0.1) * mu_reg1)
v_quiet_reg1 = np.nansum(dop.v_corr * con.image * mask.is_quiet_sun() * mu_reg1) / np.nansum(con.image * mask.is_quiet_sun() * mu_reg1)
v_umb_reg1 = np.nansum(dop.v_corr * con.image * mask.is_umbra() * mu_reg1) / np.nansum(con.image * mask.is_umbra() * mu_reg1)
v_pen_reg1 = np.nansum(dop.v_corr * con.image * mask.is_penumbra() * mu_reg1) / np.nansum(con.image * mask.is_penumbra() * mu_reg1)
v_network_reg1 = np.nansum(dop.v_corr * con.image * mask.is_network() * mu_reg1) / np.nansum(con.image * mask.is_network() * mu_reg1)
v_plage_reg1 = np.nansum(dop.v_corr * con.image * mask.is_plage() * mu_reg1) / np.nansum(con.image * mask.is_plage() * mu_reg1)


v_hat_reg2 = np.nansum(dop.v_corr * con.image * (con.mu >= 0.1) * mu_reg2) / np.nansum(con.image * (con.mu >= 0.1) * mu_reg2)
v_quiet_reg2 = np.nansum(dop.v_corr * con.image * mask.is_quiet_sun() * mu_reg2) / np.nansum(con.image * mask.is_quiet_sun() * mu_reg2)
v_umb_reg2 = np.nansum(dop.v_corr * con.image * mask.is_umbra() * mu_reg2) / np.nansum(con.image * mask.is_umbra() * mu_reg2)
v_pen_reg2 = np.nansum(dop.v_corr * con.image * mask.is_penumbra() * mu_reg2) / np.nansum(con.image * mask.is_penumbra() * mu_reg2)
v_network_reg2 = np.nansum(dop.v_corr * con.image * mask.is_network() * mu_reg2) / np.nansum(con.image * mask.is_network() * mu_reg2)
v_plage_reg2 = np.nansum(dop.v_corr * con.image * mask.is_plage() * mu_reg2) / np.nansum(con.image * mask.is_plage() * mu_reg2)


# get light in regions
f_light_reg1 = np.nansum(con.image * (con.mu >= 0.1) * mu_reg1) / all_light
f_light_reg2 = np.nansum(con.image * (con.mu >= 0.1) * mu_reg2) / all_light

pdb.set_trace()

umb_light_reg1 = np.nansum(mask.is_umbra() * con.image * mu_reg1)/all_light
pen_light_reg1 = np.nansum(mask.is_penumbra() * con.image * mu_reg1)/all_light
quiet_light_reg1 = np.nansum(mask.is_quiet_sun() * con.image * mu_reg1)/all_light
network_light_reg1 = np.nansum(mask.is_network() * con.image * mu_reg1)/all_light
plage_light_reg1 = np.nansum(mask.is_plage() * con.image * mu_reg1)/all_light

umb_light_reg2 = np.nansum(mask.is_umbra() * con.image * mu_reg2)/all_light
pen_light_reg2 = np.nansum(mask.is_penumbra() * con.image * mu_reg2)/all_light
quiet_light_reg2 = np.nansum(mask.is_quiet_sun() * con.image * mu_reg2)/all_light
network_light_reg2 = np.nansum(mask.is_network() * con.image * mu_reg2)/all_light
plage_light_reg2 = np.nansum(mask.is_plage() * con.image * mu_reg2)/all_light

f_quiet_light_reg1 = np.nansum(mask.is_quiet_sun() * con.image * mu_reg1) / np.nansum(con.image * mask.is_quiet_sun())
f_quiet_light_reg2 = np.nansum(mask.is_quiet_sun() * con.image * mu_reg2) / np.nansum(con.image * mask.is_quiet_sun())

# check the sum
f_light_reg1 == umb_light_reg1 + pen_light_reg1 + quiet_light_reg1 + network_light_reg1 + plage_light_reg1
f_light_reg2 == umb_light_reg2 + pen_light_reg2 + quiet_light_reg2 + network_light_reg2 + plage_light_reg2

# combine v_hat from regions
v_hat_combined = v_hat_reg1 * f_light_reg1 + v_hat_reg2 * f_light_reg2
v_hat_whole - v_hat_combined

# combine v_hat from regions
v_hat_combined2 = 0.0
v_hat_combined2 += v_umb_reg1 * umb_light_reg1 + v_pen_reg1 * pen_light_reg1 + v_network_reg1 * network_light_reg1 + v_plage_reg1 * plage_light_reg1 + v_quiet_reg1 * quiet_light_reg1
v_hat_combined2 += v_umb_reg2 * umb_light_reg2 + v_pen_reg2 * pen_light_reg2 + v_network_reg2 * network_light_reg2 + v_plage_reg2 * plage_light_reg2 + v_quiet_reg2 * quiet_light_reg2
v_hat_whole - v_hat_combined2

# check the v_quiet combination
v_quiet_whole1 = v_quiet_reg1 * f_quiet_light_reg1 + v_quiet_reg2 * f_quiet_light_reg2
v_quiet_whole - v_quiet_whole1

# get v_convs
v_conv_reg1 = v_hat_reg1 - v_quiet_reg1 * f_quiet_light_reg1
v_conv_umb_reg1 = v_umb_reg1 - v_quiet_reg1 * f_quiet_light_reg1
v_conv_pen_reg1 = v_pen_reg1 - v_quiet_reg1 * f_quiet_light_reg1
v_conv_network_reg1 = v_network_reg1 - v_quiet_reg1 * f_quiet_light_reg1
v_conv_plage_reg1 = v_plage_reg1 - v_quiet_reg1 * f_quiet_light_reg1

v_conv_reg2 = v_hat_reg2 - v_quiet_reg2 * f_quiet_light_reg2
v_conv_umb_reg2 = v_umb_reg2 - v_quiet_reg2 * f_quiet_light_reg2
v_conv_pen_reg2 = v_pen_reg2 - v_quiet_reg2 * f_quiet_light_reg2
v_conv_network_reg2 = v_network_reg2 - v_quiet_reg2 * f_quiet_light_reg2
v_conv_plage_reg2 = v_plage_reg2 - v_quiet_reg2 * f_quiet_light_reg2

derp = 0.0
derp += v_conv_umb_reg1 * (umb_light_reg1/f_light_reg1)
derp += v_conv_pen_reg1 * (pen_light_reg1/f_light_reg1)
derp += v_conv_network_reg1 * (network_light_reg1/f_light_reg1)
derp += v_conv_plage_reg1 * (plage_light_reg1/f_light_reg1)
derp - v_conv_reg1

derp = 0.0
derp += v_conv_umb_reg2 * (umb_light_reg2/f_light_reg2)
derp += v_conv_pen_reg2 * (pen_light_reg2/f_light_reg2)
derp += v_conv_network_reg2 * (network_light_reg2/f_light_reg2)
derp += v_conv_plage_reg2 * (plage_light_reg2/f_light_reg2)
derp - v_conv_reg2

derp = 0.0
derp += v_conv_reg1 *
derp += v_conv_reg2 *
v_conv_whole - derp

# combine v_conv from regions
v_conv_whole2 = 0.0
v_conv_whole2 += v_conv_umb_reg1 * umb_light_reg1 + v_conv_pen_reg1 * pen_light_reg1 + v_conv_network_reg1 * network_light_reg1 + v_conv_plage_reg1 * plage_light_reg1
v_conv_whole2 += v_conv_umb_reg2 * umb_light_reg2 + v_conv_pen_reg2 * pen_light_reg2 + v_conv_network_reg2 * network_light_reg2 + v_conv_plage_reg2 * plage_light_reg2


v_quiet2 = (v_quiet * umb_light + v_quiet * pen_light + v_quiet * network_light + v_quiet * plage_light)/(umb_light + pen_light + network_light + plage_light)

v_conv_umb = umb_light * (v_umb - v_quiet)
v_conv_pen = pen_light * (v_pen - v_quiet)
v_conv_network = network_light * (v_network - v_quiet)
v_conv_plage = plage_light * (v_plage - v_quiet)

v_conv_umb + v_conv_pen + v_conv_network + v_conv_plage

pdb.set_trace()
