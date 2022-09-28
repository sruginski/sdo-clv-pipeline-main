import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import os, sys, pdb, csv, glob

from sdo_pypline.paths import root
from sdo_pypline.sdo_process import *
from download_plot_data import download_plot_data

# sort out paths
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

# use style
plt.style.use(str(root) + "/" + "my.mplstyle"); plt.ioff()

def main():
    # get the sdo data to plot
    con, mag, dop, aia = download_plot_data()

    # reduce the data
    print(">>> Processing and plotting data...")
    con, mag, dop, aia, mask = reduce_sdo_images(con, mag, dop, aia, mu_thresh=0.1)

    # get intensities to plot
    aia_flat = aia.iflat[~np.isnan(aia.iflat)] / aia.ld_coeffs[0]
    con_flat = con.iflat[~np.isnan(con.iflat)] / con.ld_coeffs[0]
    mag_img = mag.image[~np.isnan(mag.image)]
    w_active = mask.w_active[~np.isnan(mag.image)]
    w_quiet = mask.w_quiet[~np.isnan(mag.image)]
    mask_dark = ~np.logical_or(mask.is_umbra(), mask.is_penumbra())[~np.isnan(mag.image)]

    # length assertion
    assert len(aia_flat) == len(con_flat) == len(mag_img)

    # get random indices
    idx1 = np.random.choice(len(aia_flat), int(np.floor(len(aia_flat)/50)))
    idx2 = np.random.choice(len(aia_flat[w_active]), int(np.floor(len(aia_flat[w_active])/50)))

    # plot elephant for AIA
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.axhline(mask.aia_thresh, ls="--", c="k")
    # ax1.axvline(24.0, ls="--", c="k")
    # ax1.scatter(np.abs(mag_img), aia_flat, s=1, alpha=0.75, color="black", rasterized=True)
    # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # ax1.set_xlabel(r"$\left| B_{r,ij} \right| \ {\rm G}$")
    # ax1.set_ylabel(r"${\rm AIA}\ 1700{\rm \AA}\ I_{{\rm flat}, ij}$")
    # fig.savefig(plotdir + "aia_flat_vs_mag.pdf", dpi=150)
    # plt.clf(); plt.close()

    # plot elephant for HMI
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.axhline(mask.con_thresh, ls="--", c="k")
    # ax1.axvline(24.0, ls="--", c="k")
    # ax1.scatter(np.abs(mag_img)[idx1], con_flat[idx1] * con.ld_coeffs[0], s=1, alpha=0.75, color="black", rasterized=True)
    # ax1.set_xscale("log")
    # # ax1.set_yscale("log")
    # ax1.set_xlabel(r"$\left| B_{ij} \right| \ {\rm G}$")
    # ax1.set_ylabel(r"${\rm HMI}\ I_{{\rm flat}, ij}$")
    # fig.savefig(plotdir + "hmi_flat_vs_mag.pdf", dpi=150)
    # plt.clf(); plt.close()

    # plot continuum vs continuum
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.axhline(mask.con_thresh/con.ld_coeffs[0], ls="--", c="k", label=r"${I_{\rm thresh, HMI}}$")
    # ax1.axvline(mask.aia_thresh/aia.ld_coeffs[0], ls="dotted", c="k", label=r"${I_{\rm thresh, AIA}}$")
    # ax1.axvline(1.57 * np.nanmean(aia_flat[w_quiet]), ls="dotted", c="k", label=r"${I_{\rm thresh, AIA}}$")
    # ax1.scatter(aia_flat[idx1], con_flat[idx1], s=1, alpha=0.75, color="black", rasterized=True, label=r"${\rm Any\ }\left| B_{r,ij}\right|$")
    # ax1.scatter(aia_flat[w_active][idx2], con_flat[w_active][idx2], s=1, alpha=0.75, color="tab:blue", rasterized=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    # # ax1.set_xscale("log")
    # # ax1.set_yscale("log")
    # ax1.set_xlabel(r"${\rm Normalized\ AIA\ 1700\ \AA\ Continuum\ Intensity}$")
    # ax1.set_ylabel(r"${\rm Normalized\ HMI\ Continuum\ Intensity}$")
    # ax1.legend()
    # fig.savefig(plotdir + "aia_vs_hmi.pdf", dpi=150)
    # plt.clf(); plt.close()

    pdb.set_trace()

    # plot the distribution of intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(mask.con_thresh1/con.ld_coeffs[0], ls="--", c="k", label=r"$0.89\ \hat{I}_{\rm quiet}$")
    ax1.axvline(mask.con_thresh2/con.ld_coeffs[0], ls=":", c="k", label=r"$0.534\ \hat{I}_{\rm quiet}$")
    ax1.axvspan(0, mask.con_thresh1/con.ld_coeffs[0], color="k", alpha=0.5)
    _, bins, patches = ax1.hist(con_flat, bins="fd", color="black", lw=2, histtype="step", density=True, label=r"${\rm Any\ }\left| B_{r}\right|$")
    _, bins, patches = ax1.hist(con_flat[w_active], bins="fd", color="tab:blue", lw=2, histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlim(6e-1, 1.25)
    ax1.set_xlabel(r"${\rm Normalized\ HMI\ Continuum\ Intensity}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12, loc="upper left")
    fig.savefig(plotdir + "fig2a.pdf")
    plt.clf(); plt.close()

    pdb.set_trace()

    # plot the distribution of dark hmi intensities
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.axvline(mask.con_thresh1/con.ld_coeffs[0], ls="--", c="k", label=r"$0.89\ \hat{I}_{\rm quiet}$")
    ax1.axvline(mask.con_thresh2/con.ld_coeffs[0], ls=":", c="k", label=r"$0.45\ \hat{I}_{\rm quiet}$")
    x1, bins, patches = ax1.hist(con_flat[con_flat < (mask.con_thresh1/con.ld_coeffs[0])], cumulative=False, bins="auto", lw=2, color="black", histtype="step", density=True)
    x2, bins, patches = ax1.hist(con_flat[(con_flat < (mask.con_thresh1/con.ld_coeffs[0])) & w_active], cumulative=False, bins="auto", lw=2, color="tab:blue", histtype="step", density=True)
    ax1.set_xscale("log")
    ax1.set_xlabel(r"${\rm HMI}\ I_{{\rm flat}, ij}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    fig.savefig(plotdir + "fig2b.pdf")
    plt.clf(); plt.close()


    # plot the distribution of intensities
    derp = np.nanmean(aia_flat[w_active & mask_dark]) - np.nanstd(aia_flat[w_active & mask_dark])
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax1.axvline(mask.aia_thresh/aia.ld_coeffs[0], ls="-.", c="k", label=r"${I_{\rm thresh, AIA}}$")
    ax1.axvline(derp, ls="-.", c="k", label=r"${I_{\rm thresh, AIA}}$")
    _, bins, patches = ax1.hist(aia_flat, bins="fd",  color="tab:purple", histtype="step", density=True, label=r"${\rm Any\ }\left| B_{r,ij}\right|$")
    _, bins, patches = ax1.hist(aia_flat[w_active & mask_dark], bins="fd", color="tab:pink", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    _, bins, patches = ax1.hist(aia_flat[w_quiet], bins="fd", color="tab:orange", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| < B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"${\rm Normalized\ AIA\ 1700\ \AA\ Continuum\ Intensity}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12)
    fig.savefig(plotdir + "continuum_dist_aia.pdf")
    plt.clf(); plt.close()

    # plot the distribution of magnetic fields
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # plt.axvline(np.percentile(np.abs(mag_img[con_flat < mask.con_thresh/con.ld_coeffs[0]]), 90), ls="dotted", c="k")
    # plt.axvline(np.percentile(np.abs(mag_img[con_flat < mask.con_thresh/con.ld_coeffs[0]]), 90), ls="dotted", c="k")
    _, bins, patches = ax1.hist(np.abs(mag_img[con_flat < mask.con_thresh/con.ld_coeffs[0]]), cumulative=False, bins="auto", color="black", histtype="step", density=True)
    _, bins, patches = ax1.hist(np.abs(mag_img[(con_flat < mask.con_thresh/con.ld_coeffs[0]) & (w_active)]), cumulative=False, bins="auto", color="tab:blue", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    _, bins, patches = ax1.hist(np.abs(mag_img[(con_flat < 0.6 * mask.con_thresh/con.ld_coeffs[0]) & (w_active)]), cumulative=False, bins="auto", color="tab:orange", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    # _, bins, patches = ax1.hist(np.abs(mag_img[(con_flat > 0.6 * mask.con_thresh/con.ld_coeffs[0]) & (con_flat < 0.89 * mask.con_thresh/con.ld_coeffs[0]) & (w_active)]), cumulative=False, bins="auto", color="tab:orange", histtype="step", density=True, label=r"$\left| B_{r,ij}\right| > B_{\rm thresh}$")
    ax1.set_xscale("log")
    ax1.set_xlabel(r"$\left| B_{r,ij} \right| \ {\rm G}$")
    ax1.set_ylabel(r"${\rm Probability\ Density}$")
    ax1.legend(fontsize=12)
    fig.savefig(plotdir + "mag_dark_dist.pdf")
    plt.clf(); plt.close()

    pdb.set_trace()


    print("derp")
    print("derp")
    print("derp")
    print("derp")
    print("derp")
    print("derp")

    # plot them
    return None

if __name__ == "__main__":
    main()


