import pdb
import numpy as np


def create_sun_mask(con, mag, dop, aia, mu_thresh=0.1):
    # calculate weights
    w_active, w_quiet = calculate_weights(mag, mu_thresh=mu_thresh)

    # set values to nan for mu less than mu_thresh
    con.mask_low_mu(mu_thresh)
    dop.mask_low_mu(mu_thresh)
    mag.mask_low_mu(mu_thresh)
    aia.mask_low_mu(mu_thresh)

    # calculate scaling factor for continuum and filtergrams
    k_hat_con = np.nansum(con.image * con.ldark * w_quiet) / np.nansum(con.ldark**2 * w_quiet)
    k_hat_aia = np.nansum(aia.image * aia.ldark * w_quiet) / np.nansum(aia.ldark**2 * w_quiet)

    # calculate velocity terms
    # TODO check dop.v_rot or dop.image - v_obs for v_phot
    # TODO check def of v_conv, subtract off v_phot or not??
    v_hat = np.nansum((dop.image - dop.v_rot - dop.v_obs) * con.image) / np.nansum(con.image)
    v_phot = np.nansum(dop.v_rot * (con.image - k_hat_con * con.ldark) * w_active) / np.nansum(con.image)
    v_quiet = np.nansum((dop.image - dop.v_rot - dop.v_obs) * con.image * w_quiet) / np.nansum(con.image * w_quiet)
    v_conv = v_hat - v_quiet - v_phot

    # calculate intensity thresholds for HMI and AIA
    con_thresh = 0.89 * np.nansum(con.iflat * w_quiet)/np.nansum(w_quiet)
    aia_thresh = np.nansum(aia.iflat * w_quiet)/np.nansum(w_quiet)

    # allocate memory for mask array
    mask = np.zeros(np.shape(con.image))

    # get thresholds for penumbrae, umbrae, and quiet sun
    ind1 = ((con.iflat <= con_thresh) & (con.iflat > (0.25 * con_thresh)))
    ind2 = (con.iflat <= (0.25 * con_thresh))
    ind3 = ((con.iflat > con_thresh))# & (aia.iflat < (1.3 * aia_thresh)))

    # set mask indices
    mask[ind1] = 1 # penumbrae
    mask[ind2] = 2 # umbrae
    mask[ind3] = 3 # quiet sun

    # bright region selection
    # ind4 = ((aia.iflat > (1.3 * aia_thresh)) & (mask != 1) & (mask != 2))
    # mask[ind4] = 4

    # make remaining regions quiet sun
    ind_rem = ((con.mu > 0.0) & (mask == 0))
    mask[ind_rem] = 3



    return v_hat, v_phot, v_quiet, v_conv
