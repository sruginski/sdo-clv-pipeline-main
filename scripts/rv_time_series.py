import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob, copy
import pandas as pd

from sdo_pypline.paths import root

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# read in the data and sort by mjd
df_vels = pd.read_csv(datadir + "processed/full_disk.csv")

# df_pixels = pd.read_csv(datadir + "pixel_stats.csv")
# df_pixels.sort_values(by=["mjd", "lo_mu"], inplace=True)
# df_pixels.drop_duplicates()

# df_light = pd.read_csv(datadir + "light_stats.csv")
# df_light.sort_values(by=["mjd", "lo_mu"], inplace=True)
# df_light.drop_duplicates()

# df_intensities = pd.read_csv(datadir + "intensities.csv")
# df_intensities.sort_values(by=["mjd"], inplace=True)
# df_intensities.drop_duplicates()

# get numbers computed for full disk
# df_vels_full = df_vels[(np.isnan(df_vels.lo_mu)) & (df_vels.region == 0.0)]
# df_pixels_full = df_pixels[np.isnan(df_pixels.lo_mu)]
# df_light_full = df_light[np.isnan(df_light.lo_mu)]


derp, idx = np.unique(np.round(df_vels.mjd), return_index=True)

# plot full disk time series
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_vels.mjd[idx], df_vels.v_hat[idx], s=2, label="v_hat")
ax1.scatter(df_vels.mjd[idx], df_vels.v_phot[idx], s=2, label="v_phot")
ax1.scatter(df_vels.mjd[idx], df_vels.v_conv[idx], s=2, label="v_conv")
ax1.scatter(df_vels.mjd[idx], df_vels.v_quiet[idx], s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity\ (m/s)}$")
ax1.legend(ncol=2, fontsize=9)
fig.savefig(plotdir + "full_disk_series.pdf")
plt.show()
plt.clf(); plt.close()

exit()

# get regions only
df_vels_regs = df_vels[df_vels.region > 0.0]

# plot quiet sun velocity by mu
df_vels_quiet_sun = df_vels_regs[df_vels_regs.region == 4]

# init fig
fig = plt.figure()
ax1 = fig.add_subplot()

# get colormarp
mus = np.unique(df_vels_quiet_sun.lo_mu)
mu_step = np.diff(mus)[0]
cbar_vals = np.linspace(mus[0] - 0.5 * mu_step, mus[-1] + 0.5 * mu_step, len(mus)+1)
cmap = copy.copy(plt.cm.rainbow)
norm = mpl.colors.BoundaryNorm(cbar_vals, cmap.N)

# loop over mu annuli
for i in mus:
    df_temp = df_vels_quiet_sun[df_vels_quiet_sun.lo_mu == i]
    ax1.scatter(df_temp.mjd, df_temp.v_quiet, color=cmap(i), s=0.5)

cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ticks=mus)
cbar.set_label(r"$\mu\ {\rm bin center}$")
cbar.ax.set_yticklabels(np.round(mus + 0.05, 2))

ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Quiet-Sun\ Velocity\ (m/s)}$")
fig.savefig(plotdir + "v_quiet_time_by_mu.pdf")
plt.clf(); plt.close()

pdb.set_trace()




# get full disk velocities
df_full = df_all[np.isnan(df_all.lo_mu)]

# get velocites for full mu annuli
df_mu = df_all[(df_all.region == 0.0) & (~np.isnan(df_all.lo_mu))]

# get velocities for regions in mu
df_regs = df_all[(df_all.region > 0.0) & (~np.isnan(df_all.lo_mu))]

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.mjd, df_full.v_hat, s=2, label="v_hat")
ax1.scatter(df_full.mjd, df_full.v_phot, s=2, label="v_phot")
ax1.scatter(df_full.mjd, df_full.v_conv, s=2, label="v_conv")
ax1.scatter(df_full.mjd, df_full.v_quiet, s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.legend(ncol=2, fontsize=9)
fig.savefig(plotdir + "full_disk_series.pdf")
plt.clf(); plt.close()

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
idx = (df_regs.region == 1) & (df_regs.v_hat != 0.0)
# ax1.scatter(df_regs.mjd[idx], df_regs.v_hat[idx], s=2, label="v_hat")
ax1.scatter(df_regs.mjd[idx], df_regs.v_phot[idx], s=2, label="v_phot")
ax1.scatter(df_regs.mjd[idx], df_regs.v_conv[idx], s=2, label="v_conv")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_quiet[idx], s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.set_title("Velocity in Spots")
ax1.legend()
fig.savefig(plotdir + "spot_series.pdf")
plt.clf(); plt.close()

# make time series
fig = plt.figure()
ax1 = fig.add_subplot()
idx = (df_regs.region == 5) & (df_regs.v_hat != 0.0)
# ax1.scatter(df_regs.mjd[idx], df_regs.v_hat[idx], s=2, label="v_hat")
ax1.scatter(df_regs.mjd[idx], df_regs.v_phot[idx], s=2, label="v_phot")
ax1.scatter(df_regs.mjd[idx], df_regs.v_conv[idx], s=2, label="v_conv")
# ax1.scatter(df_regs.mjd[idx], df_regs.v_quiet[idx], s=2, label="v_quiet")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Velocity (m/s)}$")
ax1.set_title("Velocity in Plage")
ax1.legend()
fig.savefig(plotdir + "plage_series.pdf")
plt.clf(); plt.close()

# time series of other stuff
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_stats.mjd, df_stats.ffactor, s=2, label="Mag. Filling Factor")
ax1.scatter(df_stats.mjd, df_stats.plage_frac, s=2, label="Plage")
ax1.scatter(df_stats.mjd, df_stats.network_frac, s=2, label="Network")
ax1.scatter(df_stats.mjd, df_stats.pen_frac, s=2, label="Penmbrae")
ax1.scatter(df_stats.mjd, df_stats.umb_frac, s=2, label="Umbrae")
ax1.set_xlabel(r"${\rm MJD}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_ylim(0.0, 0.12)
ax1.legend(fontsize=9)
fig.savefig(plotdir + "frac_series.pdf")
plt.clf(); plt.close()

# correlation
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_stats.ffactor, df_stats.plage_frac, s=2, label="Plage")
ax1.scatter(df_stats.ffactor, df_stats.plage_frac + df_stats.network_frac, s=2, label="Plage + Network")
ax1.set_xlabel(r"${\rm Magnetic\ Filling\ Factor}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_xlim(0.05, 0.125)
ax1.set_ylim(0.00, 0.125)
ax1.legend(fontsize=9)
fig.savefig(plotdir + "frac_corr.pdf")
plt.clf(); plt.close()
