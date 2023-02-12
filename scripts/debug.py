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
from sdo_pypline.sdo_vels import *
from download_plot_data import download_plot_data

# sort out paths
# datadir = str(root / "data") + "/"
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

def main():
    # get the directory
    datadir = "/Users/michael/Desktop/sdo_debug/"

    # get the sdo data
    con1 = datadir + "hmi_ic_45s_2014_01_03_20_01_30_tai_continuum.fits"
    mag1 = datadir + "hmi_m_45s_2014_01_03_20_01_30_tai_magnetogram.fits"
    dop1 = datadir + "hmi_v_45s_2014_01_03_20_01_30_tai_dopplergram.fits"
    aia1 = datadir + "aia_lev1_1700a_2014_01_03t20_00_30_72z_image_lev1.fits"

    con2 = datadir + "hmi_ic_45s_2014_01_07_16_01_30_tai_continuum.fits"
    mag2 = datadir + "hmi_m_45s_2014_01_07_16_01_30_tai_magnetogram.fits"
    dop2 = datadir + "hmi_v_45s_2014_01_07_16_01_30_tai_dopplergram.fits"
    aia2 = datadir + "aia_lev1_1700a_2014_01_07t16_00_30_71z_image_lev1.fits"

    con3 = datadir + "hmi_ic_45s_2014_01_11_12_01_30_tai_continuum.fits"
    mag3 = datadir + "hmi_m_45s_2014_01_11_12_01_30_tai_magnetogram.fits"
    dop3 = datadir + "hmi_v_45s_2014_01_11_12_01_30_tai_dopplergram.fits"
    aia3 = datadir + "aia_lev1_1700a_2014_01_11t12_00_30_72z_image_lev1.fits"

    # process the data
    con1, mag1, dop1, aia1, mask1 = reduce_sdo_images(con1, mag1, dop1, aia1)
    con2, mag2, dop2, aia2, mask2 = reduce_sdo_images(con2, mag2, dop2, aia2)
    con3, mag3, dop3, aia3, mask3 = reduce_sdo_images(con3, mag3, dop3, aia3)

    # get quietsun velocity
    v_quiet1 = np.nansum(dop1.v_corr * con1.image * mask1.is_quiet_sun())
    v_quiet1 /= np.nansum(con1.image * mask1.is_quiet_sun())

    # calculate the velocity in penumbrae
    region_mask1 = calc_region_mask(mask1, region=2.5)
    v_hat_whole1 = calc_velocities(con1, mag1, dop1, aia1, mask1)[0]
    v_hat1, v_phot1, _, v_conv1 = calc_velocities(con1, mag1, dop1, aia1, mask1, region_mask=region_mask1, v_quiet=v_quiet1)
    light1 = np.nansum(region_mask1 * con1.image)/np.nansum(con1.image)

    # get quietsun velocity
    v_quiet2 = np.nansum(dop2.v_corr * con2.image * mask2.is_quiet_sun())
    v_quiet2 /= np.nansum(con2.image * mask2.is_quiet_sun())

    # calculate the velocity in penumbrae
    region_mask2 = calc_region_mask(mask2, region=2.5)
    v_hat_whole2 = calc_velocities(con2, mag2, dop2, aia2, mask2)[0]
    v_hat2, v_phot2, _, v_conv2 = calc_velocities(con2, mag2, dop2, aia2, mask2, region_mask=region_mask2, v_quiet=v_quiet2)
    light2 = np.nansum(region_mask2 * con2.image)/np.nansum(con2.image)

    # get quietsun velocity
    v_quiet3 = np.nansum(dop3.v_corr * con3.image * mask3.is_quiet_sun())
    v_quiet3 /= np.nansum(con3.image * mask3.is_quiet_sun())

    # calculate the velocity in penumbrae
    region_mask3 = calc_region_mask(mask3, region=2.5)
    v_hat_whole3 = calc_velocities(con3, mag3, dop3, aia3, mask3)[0]
    v_hat3, v_phot3, _, v_conv3 = calc_velocities(con3, mag3, dop3, aia3, mask3, region_mask=region_mask3, v_quiet=v_quiet3)
    light3 = np.nansum(region_mask2 * con3.image)/np.nansum(con3.image)



    cmap = plt.get_cmap("seismic").copy()
    cmap.set_bad(color="none")


    # test the doppler correction
    derp5 = dop2.image - dop2.v_obs - dop2.v_grav - dop1.v_rot - dop1.v_mer - dop2.v_cbs
    plt.imshow(derp5, cmap=cmap, vmin=-2500, vmax=2500)
    plt.colorbar()
    plt.show()

    # do some plotting
    derp0 = np.copy(dop2.v_corr)
    derp1 = np.copy(dop1.v_corr)
    derp2 = np.copy(dop2.v_corr)
    derp3 = np.copy(dop3.v_corr)

    derp4 = derp2 - derp1

    derp1[~mask1.is_penumbra()] = np.nan
    derp2[~mask2.is_penumbra()] = np.nan
    derp3[~mask3.is_penumbra()] = np.nan
    derp4[~np.logical_or(mask1.is_penumbra(), mask2.is_penumbra())] = np.nan

    plt.imshow(con1.mu, cmap="Greys")
    plt.imshow(derp1, cmap=cmap, vmin=-2500, vmax=2500)
    plt.imshow(derp2, cmap=cmap, vmin=-2500, vmax=2500)
    plt.imshow(derp3, cmap=cmap, vmin=-2500, vmax=2500)
    plt.grid(False)
    plt.colorbar()
    plt.show()

    plt.imshow(con1.mu, cmap="Greys")
    plt.imshow(derp2 - np.nanmean(derp2), cmap=cmap, vmin=-2500, vmax=2500)
    plt.grid(False)
    plt.colorbar()
    plt.show()

    plt.imshow(con1.mu, cmap="Greys")
    plt.imshow(derp4, cmap=cmap, vmin=-2500, vmax=2500)
    plt.grid(False)
    plt.colorbar()
    plt.show()


    pdb.set_trace()



if __name__ == "__main__":
    main()
