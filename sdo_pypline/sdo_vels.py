import numpy as np
import pdb

def calc_region_mask(mask, region=None, hi_mu=None, lo_mu=None):
    # get mask for region type specified
    if region is None:
        region_mask = True
    elif region in np.unique(mask.regions):
        region_mask = (region == mask.regions)

    # get masks for mu annuli
    if ((hi_mu is None) | (lo_mu is None)):
        region_mask *= True
    else:
        assert lo_mu < hi_mu
        region_mask *= ((mask.mu > lo_mu) & (mask.mu <= hi_mu))

    return region_mask

def calc_velocities(con, mag, dop, aia, mask, region_mask=True, weight_denom=True):
    # don't bother doing math if there is nothing in the mask
    if (type(region_mask) is np.ndarray) and (~region_mask.any()):
        return 0.0, 0.0, 0.0, 0.0

    # get weights
    w_quiet = mask.w_quiet
    w_active = mask.w_active

    # calculate scaling factor for continuum and filtergrams
    k_hat_con = np.nansum(con.image * con.ldark * w_quiet) / np.nansum(con.ldark**2 * w_quiet)

    # calculate velocity terms
    v_hat = np.nansum((dop.image - dop.v_rot - dop.v_obs) * con.image * region_mask)
    v_phot = np.nansum(dop.v_rot * (con.image - k_hat_con * con.ldark) * w_active * region_mask)

    # get quiet sun velocity
    if (np.nansum(w_quiet * region_mask) == 0):
        v_quiet = 0.0
    else:
        v_quiet = np.nansum((dop.image - dop.v_rot - dop.v_obs) * con.image * w_quiet * region_mask)

    # divide velocities by the denominator
    if weight_denom:
        v_hat /= np.nansum(con.image * region_mask)
        v_phot /= np.nansum(con.image * region_mask)
        if v_quiet != 0.0:
            v_quiet /= np.nansum(con.image * w_quiet * region_mask)
    else:
        v_hat /= np.nansum(con.image)
        v_phot /= np.nansum(con.image)
        if v_quiet != 0.0:
            v_quiet /= np.nansum(con.image * w_quiet)

    # get convective velocity by subtracting off other terms
    v_conv = v_hat - v_quiet

    return v_hat, v_phot, v_quiet, v_conv

def calc_mag_stats(con, mag, mask, region_mask=True, weight_denom=True):
    # don't bother doing math if there is nothing in the mask
    if (type(region_mask) is np.ndarray) and (~region_mask.any()):
        return 0.0, 0.0, 0.0

    # get average and std for mag field strength
    abs_mag = np.abs(mag.image) * region_mask
    mag_avg = np.nanmean(abs_mag)
    mag_std = np.nanstd(abs_mag)

    # get intensity weighted unsigned magnetic field strength
    mag_unsigned = np.nansum(np.abs(mag.image) * con.image * region_mask)

    # divide by the denominator
    if weight_denom:
        mag_unsigned /= np.nansum(con.image * region_mask)
    else:
        mag_unsigned /= np.nansum(con.image)

    return mag_avg, mag_std, mag_unsigned
