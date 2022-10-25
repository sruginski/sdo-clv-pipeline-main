import numpy as np
import pdb, warnings
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from sunpy.map import Map as sun_map
from sunpy.coordinates import frames

from scipy import ndimage
from astropy.wcs import WCS
from scipy.optimize import curve_fit
from reproject import reproject_interp
from skimage.measure import regionprops
from astropy.wcs import FITSFixedWarning
from astropy.io.fits.verify import VerifyWarning

from .sdo_io import *
from .limbdark import *
from .legendre import *

warnings.simplefilter("ignore", category=VerifyWarning)
warnings.simplefilter("ignore", category=FITSFixedWarning)

class SDOImage(object):
    def __init__(self, file):
        # set the filename
        self.filename = file

        # set attributes from the heaader
        self.parse_header(file)

        # read in the data
        self.image = read_data(file)

        # initialize mu_thresh
        self.mu_thresh = 0.0
        return None

    def parse_header(self, file):
        # read the header
        head = read_header(file)
        self.wcs = WCS(head)

        # parse it
        self.naxis1 = head["NAXIS1"]
        self.naxis2 = head["NAXIS2"]
        self.wavelength = head["WAVELNTH"]
        self.crpix1 = head["CRPIX1"]
        self.crpix2 = head["CRPIX2"]
        self.crval1 = head["CRVAL1"]
        self.crval2 = head["CRVAL2"]
        self.cdelt1 = head["CDELT1"]
        self.cdelt2 = head["CDELT2"]
        self.crota2 = np.deg2rad(head["CROTA2"])
        self.dsun_obs = head["DSUN_OBS"]
        self.dsun_ref = head["DSUN_REF"]
        self.rsun_obs = head["RSUN_OBS"]
        self.rsun_ref = head["RSUN_REF"]
        self.crln_obs = np.deg2rad(head["CRLN_OBS"])
        self.crlt_obs = np.deg2rad(head["CRLT_OBS"])
        self.car_rot = head["CAR_ROT"]
        self.obs_vr = head["OBS_VR"]
        self.obs_vw = head["OBS_VW"]
        self.obs_vn = head["OBS_VN"]
        self.rsun_obs = head["RSUN_OBS"]
        self.date_obs = head["DATE-OBS"]
        if "CONTENT" in head.keys():
            self.content = head["CONTENT"]
        else:
            self.content = "FILTERGRAM"
        return None

    def calc_geometry(self):
        # get distance to sun in solar radii and focal length
        self.dist_sun = self.dsun_obs/self.rsun_ref
        self.focal_len = 180. * 3600. / np.pi / self.cdelt1

        # mesh of pixels and distances to pixels in pixel units
        paxis1 = np.linspace(1, self.naxis1, self.naxis1) - self.crpix1
        paxis2 = np.linspace(1, self.naxis2, self.naxis2) - self.crpix2
        self.px, self.py = np.meshgrid(paxis1, paxis2, sparse=True)
        self.pr = np.sqrt(self.px**2.0 + self.py**2.0)

        #  distances in solar radii
        rr_complex = self.focal_len**2 * self.pr**2 + self.pr**4 - self.dist_sun**2 * self.pr**4 + 0j
        self.rr = np.real(self.dist_sun * self.focal_len * self.pr - np.sqrt(rr_complex)) / (self.focal_len**2 + self.pr**2)
        self.rr_obs = np.real(np.sqrt(1.0 - self.rr**2 + 0j))

        # calculate grid of mus
        cos_alpha = self.rr_obs
        sin_alpha = self.rr
        cos_theta = (self.dist_sun - cos_alpha) / np.sqrt(self.rr**2 + (self.dist_sun - cos_alpha)**2)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        self.mu = np.real(cos_alpha * cos_theta - sin_alpha * sin_theta)
        return None

    def inherit_geometry(self, other_image):
        self.dist_sun = np.copy(other_image.dist_sun)
        self.focal_len = np.copy(other_image.focal_len)
        self.px = np.copy(other_image.px)
        self.py = np.copy(other_image.py)
        self.pr = np.copy(other_image.pr)
        self.rr = np.copy(other_image.rr)
        self.rr_obs = np.copy(other_image.rr_obs)
        self.mu = np.copy(other_image.mu)
        return None

    def is_magnetogram(self):
        return self.content == "MAGNETOGRAM"

    def is_dopplergram(self):
        return self.content == "DOPPLERGRAM"

    def is_continuum(self):
        return self.content == "CONTINUUM INTENSITY"

    def is_filtergram(self):
        return self.content == "FILTERGRAM"

    def mask_low_mu(self, mu_thresh):
        self.mu_thresh = mu_thresh
        self.image[np.logical_or(self.mu <= mu_thresh, np.isnan(self.mu))] = np.nan

        if self.is_continuum() | self.is_filtergram():
            self.ldark[np.logical_or(self.mu <= mu_thresh, np.isnan(self.mu))] = np.nan
            self.iflat[np.logical_or(self.mu <= mu_thresh, np.isnan(self.mu))] = np.nan
        return None

    def correct_magnetogram(self):
        assert self.is_magnetogram()
        self.image /= self.mu
        return None

    def correct_dopplergram(self):
        assert self.is_dopplergram()

        # read file in as sunpy map to do coordinate transforms
        # methods adapted from https://arxiv.org/abs/2105.12055
        # original implementation at https://github.com/samarth-kashyap/hmi-clean-ls
        smap = sun_map(self.filename)

        # get coordinates
        paxis1 = np.arange(self.naxis1)
        paxis2 = np.arange(self.naxis2)
        xx, yy = np.meshgrid(paxis1, paxis2) * u.pix
        hpc = smap.pixel_to_world(xx, yy)    # helioprojective cartesian
        hgs = hpc.transform_to(frames.HeliographicStonyhurst)
        hcc = hpc.transform_to(frames.Heliocentric) # heliocentric cartiesian

        # get coordinates
        self.xx = hcc.x
        self.yy = hcc.y
        self.rr = np.sqrt(hpc.Tx**2 + hpc.Ty**2) / smap.rsun_obs
        self.B0 = smap.observer_coordinate.lat
        self.lat = hgs.lat + 90 * u.deg
        self.lon = hgs.lon

        # get distance to sun in solar radii
        self.rsun_solrad = smap.observer_coordinate.radius/smap.rsun_meters

        # get mask excluding nans / sqrts of negatives
        self.rr[self.rr > 0.95] = np.nan
        self.mask_nan = np.logical_and(~np.isnan(self.rr), ~np.isnan(self.lat))

        # velocity components
        self.v_grav = 632 # m/s, constant
        self.calc_spacecraft_vel() # spacecraft velocity
        self.calc_bulk_vel() # differential rotation + meridional flows + cbs
        return None

    def calc_spacecraft_vel(self):
        assert self.is_dopplergram()

        # pre-compute trigonometric quantities
        sig = np.arctan(self.rr/self.rsun_solrad)
        chi = np.arctan2(self.xx, self.yy)
        sin_sig = np.sin(sig)
        cos_sig = np.cos(sig)
        sin_chi = np.sin(chi)
        cos_chi = np.cos(chi)

        # get satellite velocities
        self.obs_vw
        self.obs_vn
        self.obs_vr

        # project satellite velocity into coordinate frame
        vr1 = self.obs_vr * cos_sig
        vr2 = -self.obs_vw * sin_sig * sin_chi
        vr3 = -self.obs_vn * sin_sig * cos_chi

        # reshape into 4096 x 4096 array
        self.v_obs = np.zeros((4096, 4096))
        self.v_obs[self.mask_nan] = -(vr1 + vr2 + vr3)[self.mask_nan]
        self.v_obs[~self.mask_nan] = np.nan
        return None

    def calc_bulk_vel(self):
        # methods adapted from https://arxiv.org/abs/2105.12055
        # original implementation at https://github.com/samarth-kashyap/hmi-clean-ls
        assert self.is_dopplergram()

        # pre-compute trigonometric quantities
        cos_B0 = np.cos(self.B0)
        sin_B0 = np.sin(self.B0)

        self.lat_mask = self.lat[self.mask_nan].copy()
        self.lon_mask = self.lon[self.mask_nan].copy()
        self.rho_mask = self.rr[self.mask_nan].copy()

        cos_theta = np.cos(self.lat_mask)
        sin_theta = np.sin(self.lat_mask)
        cos_phi = np.cos(self.lon_mask)
        sin_phi = np.sin(self.lon_mask)

        self.lt = sin_B0 * sin_theta - cos_B0 * cos_theta * cos_phi
        self.lp = cos_B0 * sin_phi

        # calculate legendre poylnomials
        print(">>> Generating ~Legendre~ Polynomials")
        pl_theta, dt_pl_theta = gen_leg(5, self.lat_mask)
        pl_rho, dt_pl_rho = gen_leg_x(5, self.rho_mask)

        # allocate memory for linalg operations
        n_poly = 11
        self.im_arr = np.zeros((n_poly, self.lt.shape[0]))
        A = np.zeros((n_poly, n_poly))

        # differential rotation (axisymmetric feature; s = 1, 3, 5)
        self.im_arr[0, :] = dt_pl_theta[1, :] * self.lp
        self.im_arr[1, :] = dt_pl_theta[3, :] * self.lp
        self.im_arr[2, :] = dt_pl_theta[5, :] * self.lp

        # meridional circulation (axisymmetric feature; s = 2, 4)
        self.im_arr[3, :] = dt_pl_theta[2, :] * self.lt
        self.im_arr[4, :] = dt_pl_theta[4, :] * self.lt

        # axisymmetric feature (frame=pole at disk-center)
        # s = 0-5
        self.im_arr[5, :] = pl_rho[0, :]
        self.im_arr[6, :] = pl_rho[1, :]
        self.im_arr[7, :] = pl_rho[2, :]
        self.im_arr[8, :] = pl_rho[3, :]
        self.im_arr[9, :] = pl_rho[4, :]
        self.im_arr[10, :] = pl_rho[5, :]

        # do the linear algebra
        self.dat = (self.image - self.v_obs - self.v_grav)[self.mask_nan].copy()
        self.RHS = self.im_arr.dot(self.dat)

        for i in range(n_poly):
            for j in range(n_poly):
                A[i, j] = self.im_arr[i, :].dot(self.im_arr[j, :])

        Ainv = inv_SVD(A, 1e5)
        self.fit_params = Ainv.dot(self.RHS)

        # get total effect
        self.v_idk = np.zeros((4096, 4096))
        self.v_idk[self.mask_nan] = self.fit_params.dot(self.im_arr)
        self.v_idk[~self.mask_nan] = np.nan

        # get rotation
        self.v_rot = np.zeros((4096, 4096))
        self.v_rot[self.mask_nan] = self.fit_params[:3].dot(self.im_arr[:3, :])
        self.v_rot[~self.mask_nan] = np.nan

        # get meridional circulation
        self.v_mer = np.zeros((4096, 4096))
        self.v_mer[self.mask_nan] = self.fit_params[3:5].dot(self.im_arr[3:5, :])
        self.v_mer[~self.mask_nan] = np.nan

        # get convective blueshift w/ limb
        self.v_cbs = np.zeros((4096, 4096))
        self.v_cbs[self.mask_nan] = self.fit_params[5:].dot(self.im_arr[5:, :])
        self.v_cbs[~self.mask_nan] = np.nan

        # get corrected velocity
        self.dat -= self.fit_params.dot(self.im_arr)
        self.v_corr = np.zeros((4096, 4096))
        self.v_corr[self.mask_nan] = self.dat
        self.v_corr[~self.mask_nan] = np.nan
        return None

    def calc_limb_darkening(self, mu_lim=0.1, num_mu=50):
        assert (self.is_continuum() | self.is_filtergram())

        # get average intensity in evenly spaced rings
        mu_edge = np.linspace(1.0, mu_lim, num=num_mu)
        avg_int = np.zeros(len(mu_edge)-1)
        for i in range(len(avg_int)):
            # find indices in ring that aren't nan
            inds = (self.mu > mu_edge[i+1]) & (self.mu <= mu_edge[i]) & (~np.isnan(self.image))

            # mask section that are big outliers
            ints = self.image[inds]
            ints[np.abs(ints - np.mean(ints)) >= (3.0 * np.std(ints))] = np.nan
            avg_int[i] = np.nanmean(ints)

        # take averages in mu annuli to fit to
        mu_avgs = (mu_edge[1:] + mu_edge[0:-1]) / 2.0

        # set the initial guess parameters for optimization
        if self.is_continuum():
            p0 = [59000.0, 0.38, 0.23]
        elif self.is_filtergram():
            p0 = [1000, 0.9, -0.25]
        else:
            return None

        # do the fit and divide out the LD profile
        popt, pcov = curve_fit(quad_darkening, mu_avgs, avg_int, p0=p0)
        self.ld_coeffs = popt
        self.ldark = quad_darkening_two(self.mu, *popt[1:])
        self.iflat = self.image/self.ldark
        return None

    def rescale_to_hmi(self, hmi_image):
        assert self.is_filtergram()

        # rescale the image
        self.image = reproject_interp((self.image, read_header(self.filename)),
                                      read_header(hmi_image.filename),
                                      return_footprint=False)

        # TODO recalculate the mu???
        self.wcs = hmi_image.wcs
        self.mu = hmi_image.mu

# for creating pixel mask with thresholded regions
def calculate_weights(mag):
    # set magnetic threshold
    mag_thresh = 24.0/mag.mu

    # make flag array for magnetically active areas
    w_active = (np.abs(mag.image) > mag_thresh).astype(float)

    # convolve with boxcar filter to remove isolated pixels
    w_conv = ndimage.convolve(w_active, np.ones([3,3]), mode="constant")
    w_active = np.logical_and(w_conv >= 2., w_active == 1.)
    w_active[np.logical_or(mag.mu <= mag.mu_thresh, np.isnan(mag.mu))] = False

    # make weights array for magnetically quiet areas
    w_quiet = ~w_active
    w_quiet[np.logical_or(mag.mu <= mag.mu_thresh, np.isnan(mag.mu))] = False
    return w_active, w_quiet

class SunMask(object):
    def __init__(self, con, mag, dop, aia):
        # check argument order/names are correct
        assert con.is_continuum()
        assert mag.is_magnetogram()
        assert dop.is_dopplergram()
        assert aia.is_filtergram()

        # copy observation date
        self.date_obs = con.date_obs

        # inherit the geometry and the WCS
        self.wcs = WCS(read_header(con.filename))
        self.inherit_geometry(con)

        # calculate weights
        self.w_active, self.w_quiet = calculate_weights(mag)

        # calculate magnetic filling factor
        npix = np.nansum(con.mu > con.mu_thresh)
        self.ff = np.nansum(self.w_active[con.mu > con.mu_thresh]) / npix

        # identify regions
        self.identify_regions(con, mag, dop, aia)

        # get region fracs
        self.pen_frac = np.nansum(self.is_penumbra()) / npix
        self.umb_frac = np.nansum(self.is_umbra()) / npix
        self.quiet_frac = np.nansum(self.is_quiet()) / npix
        self.network_frac = np.nansum(self.is_network()) / npix
        self.plage_frac = np.nansum(self.is_plage()) / npix

        return None

    def inherit_geometry(self, other_image):
        self.dist_sun = other_image.dist_sun
        self.focal_len = other_image.focal_len
        self.px = other_image. px
        self.py = other_image.py
        self.pr = other_image.pr
        self.rr = other_image.rr
        self.rr_obs = other_image.rr_obs
        self.mu = other_image.mu
        return None

    def identify_regions(self, con, mag, dop, aia):
        # allocate memory for mask array
        self.regions = np.zeros(np.shape(con.image))

        # calculate intensity thresholds for HMI
        self.con_thresh1 = 0.89 * np.nansum(con.iflat * self.w_quiet)/np.nansum(self.w_quiet)
        self.con_thresh2 = 0.45 * np.nansum(con.iflat * self.w_quiet)/np.nansum(self.w_quiet)

        # get indices for umbrae
        ind1 = con.iflat <= self.con_thresh2

        # get indices for penumbrae (blueshifted and redshifted)
        indp = (con.iflat <= self.con_thresh1) & (con.iflat > self.con_thresh2)
        ind2 = (dop[indp] <= 0.0)
        ind3 = (dop[indp] > 0.0)

        # get indices for quiet sun
        ind4 = (con.iflat > self.con_thresh1) & self.w_quiet

        # calculate intensity thresholds for AIA
        weights = self.w_active * (~ind1) * (~ind2) * (~ind3)
        self.aia_thresh = np.nansum(aia.iflat * weights)/np.nansum(weights)

        # get indices for bright regions (plage/faculae + network)
        ind5a = (con.iflat > self.con_thresh1) & self.w_active
        ind5b = (aia.iflat > self.aia_thresh) & (~ind1) & (~ind2) & (~ind3)
        ind5 = ind5a | ind5b

        # set mask indices
        self.regions[ind1] = 1 # umbrae
        self.regions[ind2] = 2 # blueshifted penumbrae
        self.regions[ind2] = 3 # redshifted penumbrae
        self.regions[ind4] = 4 # quiet sun
        self.regions[ind5] = 5 # bright areas (will separate into plage + network)

        # label unique contiguous bright regions and calculate their sizes
        binary_img = self.regions == 5
        structure = ndimage.generate_binary_structure(2,2)
        labels, nlabels = ndimage.label(binary_img, structure=structure)

        # get labeled region areas and perimeters
        rprops = regionprops(labels)
        areas = np.array([rprop.area for rprop in rprops]).astype(float)
        areas *= (1e6/np.sum(self.mu > 0.0)) # convert to microhemispheres

        # calculate the perimeter to area ratio and apply threshold
        # TODO: need to handle projection/foreshortening, or is pixel count fine
        area_thresh = 20

        # assign region type to plage for ratios less than ratio thresh
        ind6 = np.concatenate(([False], areas > area_thresh))[labels]
        self.regions[ind6] = 6 # plage

        # set isolated bright pixels to quiet sun
        ind_iso = np.concatenate(([False], areas == 1))[labels]
        self.regions[ind_iso] = 4 # quiet sun

        # make any remaining unclassified pixels quiet sun
        ind_rem = ((con.mu > con.mu_thresh) & (self.regions == 0))
        self.regions[ind_rem] = 4 # quiet sun

        # set values beyond mu_thresh to nan
        self.regions[np.logical_or(con.mu <= con.mu_thresh, np.isnan(con.mu))] = np.nan

        return None

    def is_umbra(self):
        return self.regions == 1

    def is_penumbra(self):
        return np.logical_or(self.regions == 2, self.regions == 3)

    def is_quiet(self):
        return self.regions == 4

    def is_network(self):
        return self.regions == 5

    def is_plage(self):
        return self.regions == 6
