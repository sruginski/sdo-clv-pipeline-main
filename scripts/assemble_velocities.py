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

def calc_whole_disk_vel(df_vels, df_light):
    # get whole disk
    v_hat_whole = df_vels.v_hat[0]
    v_phot_whole = df_vels.v_phot[0]
    v_quiet_whole = df_vels.v_quiet[0]
    v_conv_whole = df_vels.v_conv[0]

    # get mu slices
    df_vels_annuli = df_vels.iloc[1:]
    df_light_annuli = df_light.iloc[1:]

    # lo_mus
    lo_mus = np.unique(df_vels_annuli.lo_mu)
    hi_mus = np.unique(df_vels_annuli.hi_mu)

    # loop over mu annuli
    v_hat_whole2 = 0.0
    v_phot_whole2 = 0.0
    v_quiet_whole2 = 0.0
    v_conv_whole2 = 0.0
    v_hat_mu = np.zeros(len(lo_mus))
    v_phot_mu = np.zeros(len(lo_mus))
    v_quiet_mu = np.zeros(len(lo_mus))
    v_conv_mu = np.zeros(len(lo_mus))
    fracs_mu = np.zeros(len(lo_mus))
    for i in range(len(lo_mus)):
        # get just this mu_annulus
        df_vels_mu = df_vels_annuli[df_vels_annuli.lo_mu == lo_mus[i]]
        df_light_mu = df_light_annuli[df_light_annuli.lo_mu == lo_mus[i]]

        fracs = [df_light_mu.umb_frac.item(), df_light_mu.blu_pen_frac.item(), df_light_mu.red_pen_frac.item(), df_light_mu.quiet_frac.item(), df_light_mu.network_frac.item(), df_light_mu.plage_frac.item()]

        # do the calculation
        for j, k in enumerate(np.unique(df_vels_mu.region)):
            # get subset
            k = int(k)
            df_vels_reg = df_vels_mu[df_vels_mu.region == k]

            # get v_quiet for mu
            v_quiet_temp = df_vels_mu[df_vels_mu.region == 4].v_quiet.item()
            v_quiet_frac = fracs[3]#/np.sum(fracs)

            # add em up
            fracs_mu[i] += fracs[j]
            v_hat_whole2 += (df_vels_reg.v_hat.item() * fracs[j])
            v_hat_mu[i] += (df_vels_reg.v_hat.item() * fracs[j])
            v_phot_whole2 += (df_vels_reg.v_phot.item() * fracs[j])
            v_phot_mu[i] += (df_vels_reg.v_phot.item() * fracs[j])
            v_quiet_whole2 += (df_vels_reg.v_quiet.item() * fracs[j])
            v_quiet_mu[i] += (df_vels_reg.v_quiet.item() * fracs[j])

            # pdb.set_trace()
            # derp = df_vels_reg.v_hat.item() - (v_quiet_temp * v_quiet_frac)
            derp = df_vels_reg.v_conv.item()

            v_conv_whole2 += derp
            v_conv_mu[i] += derp

    print(v_conv_whole)
    print(v_conv_whole2)
    pdb.set_trace()


    print("derp")

    return None

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# read in the data and sort by mjd
df_all = pd.read_csv(datadir + "velocities.csv")
df_all.sort_values(by=["mjd"], inplace=True)

df_pixels = pd.read_csv(datadir + "pixel_stats.csv")
df_pixels.sort_values(by=["mjd"], inplace=True)

df_light = pd.read_csv(datadir + "light_stats.csv")
df_light.sort_values(by=["mjd"], inplace=True)

for mjd in np.unique(df_all.mjd):
    row_vels = df_all[df_all.mjd == mjd]
    row_light = df_light[df_light.mjd == mjd]

    calc_whole_disk_vel(row_vels, row_light)

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
