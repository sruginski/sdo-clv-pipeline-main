import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob
import pandas as pd

# use style
plt.style.use("my.mplstyle"); plt.ioff()

# data files
# datdir = "/Users/michael/Desktop/sdo-pypline/data/"
datdir = "/Users/michael/Desktop/sdo_output/"
df_full = pd.read_csv(datdir + "rv_full_disk.csv")
df_regs = pd.read_csv(datdir + "rv_regions.csv")
df_mu = pd.read_csv(datdir + "rv_mu.csv")
df_aia_ld = pd.read_csv(datdir + "aia_ld_params.csv")
df_hmi_ld = pd.read_csv(datdir + "hmi_ld_params.csv")
df_con_thresh = pd.read_csv(datdir + "con_thresh.csv")
df_mag_stats = pd.read_csv(datdir + "mag_stats.csv")

# correlation
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(df_full.ffactor, df_full.plage_frac, s=2, label="Plage")
ax1.scatter(df_full.ffactor, df_full.plage_frac + df_full.network_frac, s=2, label="Plage + Network")
ax1.set_xlabel(r"${\rm Magnetic\ Filling\ Factor}$")
ax1.set_ylabel(r"${\rm Fraction}$")
ax1.set_xlim(0.05, 0.125)
ax1.set_ylim(0.00, 0.125)
ax1.legend(fontsize=9)
fig.savefig(datdir + "frac_corr.pdf")
plt.clf(); plt.close()
