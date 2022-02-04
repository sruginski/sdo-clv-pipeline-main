def create_sun_mask(con, mag, dop, aia, mu_thresh=0.1):
    # set magnetic threshold
    mag_thresh = 24.0/mag.mu

    # make flag array for magnetically active areas
    weights_active = np.zeros(np.shape(mag.image))
    weights_active[np.abs(mag.image) > mag_thresh] = 1.0
    weights_active[mag.mu <= mu_thresh] = 0.0

    # make weights array for magnetically quiet areas
    weights_quiet = np.copy(weights_active)
    weights_quiet[weights_active == 1.0] = 0.0
    weights_quiet[weights_active == 0.0] = 1.0
    weights_quiet[mag.mu <= mu_thresh] = 0.0

    # calculate intensity threshold
    on_mu = con.mu > mu_thresh
    con_thresh = 0.89 * sum(con.image[on_mu] * weights_quiet[on_mu])/sum(weights_quiet[on_mu])

    # calculate AIA threshold
    on_mu = aia.mu > mu_thresh
    aia_thresh = sum(aia.image[on_mu] * weights_quiet[on_mu])/sum(weights_quiet[on_mu])


    return
