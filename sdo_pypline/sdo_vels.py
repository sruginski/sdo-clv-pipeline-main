import numpy as np
import pdb

def calc_velocities(con, mag, dop, aia, mask, region=None, hi_mu=None, lo_mu=None):
    # get weights
    w_quiet = mask.w_quiet
    w_active = mask.w_active

    # get mask for region type specified
    if region is None:
        region_mask = True
    elif region in np.unique(mask.regions):
        region_mask = (region == mask.regions)

    # get masks for mu annuli
    if ((hi_mu is None) | (lo_mu is None)):
        mu_mask = True
    else:
        assert lo_mu < hi_mu
        mu_mask = ((con.mu > lo_mu) & (con.mu <= hi_mu))


    # calculate scaling factor for continuum and filtergrams
    k_hat_con = np.nansum(con.image * con.ldark * w_quiet) / np.nansum(con.ldark**2 * w_quiet)
    # k_hat_aia = np.nansum(aia.image * aia.ldark * w_quiet) / np.nansum(aia.ldark**2 * w_quiet)

    # calculate velocity terms
    # TODO check dop.v_rot or dop.image - v_obs for v_phot
    # TODO check def of v_conv, subtract off v_phot or not??
    v_hat = np.nansum((dop.image - dop.v_rot - dop.v_obs) * con.image * region_mask) / np.nansum(con.image * region_mask)
    v_phot = np.nansum(dop.v_rot * (con.image - k_hat_con * con.ldark) * w_active * region_mask) / np.nansum(con.image * region_mask)

    # get quiet sun velocity
    if (np.nansum(w_quiet * region_mask) == 0):
        v_quiet = 0.0
    else:
        v_quiet = np.nansum((dop.image - dop.v_rot - dop.v_obs) * con.image * w_quiet * region_mask) / np.nansum(con.image * w_quiet * region_mask)

    # get convective velocity by subtracting off other terms
    v_conv = v_hat - v_quiet - v_phot

    return v_hat, v_phot, v_quiet, v_conv
