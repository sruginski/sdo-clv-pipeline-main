import numpy as np
import pdb

def calc_region_mask(mask, region=None, hi_mu=None, lo_mu=None):
    # get mask for region type specified
    if region is None:
        region_mask = True
    elif region in np.unique(mask.regions):
        region_mask = ((region == mask.regions) & (mask.mu >= mask.mu_thresh))

    # get masks for mu annuli
    if ((hi_mu is None) | (lo_mu is None)):
        region_mask *= True
    else:
        assert lo_mu < hi_mu
        region_mask *= ((mask.mu > lo_mu) & (mask.mu <= hi_mu)  & (mask.mu >= mask.mu_thresh))

    return region_mask

def calc_velocities(con, mag, dop, aia, mask, region_mask=None, v_quiet=None):
    # don't bother doing math if there is nothing in the mask
    if (type(region_mask) is np.ndarray) and (~region_mask.any()):
        return 0.0, 0.0, 0.0, 0.0

    # get default region_mask
    if region_mask is None:
        region_mask = (con.mu >= con.mu_thresh)

    # get weights
    w_quiet = mask.is_quiet()
    w_active = ~w_quiet

    # calculate scaling factor for continuum and filtergrams
    k_hat_con = np.nansum(con.image * con.ldark * w_quiet) / np.nansum(con.ldark**2 * w_quiet)

    # calculate velocity terms
    v_hat = np.nansum(dop.v_corr * con.image * region_mask)

    # TODO add v_mer???
    v_phot = np.nansum(dop.v_rot * (con.image - k_hat_con * con.ldark) * w_active * region_mask)

    # divide velocities by the denominator (only calculate it once)
    denom = np.nansum(con.image * region_mask)
    v_hat /= denom
    v_phot /= denom

    if v_quiet is None:
        # calculate v_quiet
        v_quiet = np.nansum(dop.v_corr * con.image * w_quiet * region_mask)

        # TODO change of definition - divide out quiet light or total light?
        # v_quiet /= np.nansum(con.image * region_mask)
        v_quiet /= np.nansum(con.image * w_quiet * region_mask)

        # get convective velocity by subtracting off other terms
        v_conv = v_hat - v_quiet # - v_phot
        return v_hat, v_phot, v_quiet, v_conv
    else:
        # get convective velocity by subtracting off other terms
        v_conv = v_hat - v_quiet # - v_phot
        return v_hat, v_phot, 0.0, v_conv
    return None

def calc_mag_stats(con, mag, mask, region_mask=True):
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
    mag_unsigned /= np.nansum(con.image * region_mask)

    return mag_avg, mag_std, mag_unsigned
