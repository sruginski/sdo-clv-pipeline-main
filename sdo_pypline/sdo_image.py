import numpy as np
from .sdo_io import *
from .limbdark import *

import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy import ndimage
from scipy.optimize import curve_fit
from reproject import reproject_interp

class SDOImage:
    def __init__(self, file):
        # set attributes from the heaader
        self.parse_header(file)

        # read in the data
        self.image = read_data(file)

        # initialize mu_thresh
        self.mu_thresh = 0.0

        # get distance to sun in solar radii and focal length
        self.dist_sun = self.dsun_obs/self.rsun_ref
        self.focal_len = 180. * 3600. / np.pi / self.cdelt1

        # mesh of pixels and distances to pixels in pixel units
        paxis1 = np.arange(self.naxis1) - self.crpix1
        paxis2 = np.arange(self.naxis2) - self.crpix2
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

    def parse_header(self, file):
        # read the header
        head = read_header(file)
        self.head = head

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

    def calc_vrot_vobs(self):
        assert self.is_dopplergram()

        # pre-compute geometric quantities
        coscrlt  = np.cos(self.crlt_obs)
        sincrlt  = np.sin(self.crlt_obs)
        coscrln  = np.cos(self.crln_obs)
        sincrln  = np.sin(self.crln_obs)
        coscrota = np.cos(self.crota2)
        sincrota = np.sin(self.crota2)

        # heliographic coordinates
        dw = self.py * sincrota + self.px * coscrota
        dn = self.py * coscrota - self.px * sincrota

        # cartesian coordinates
        rw_obs = self.rr * dw/self.pr
        rn_obs = self.rr * dn/self.pr
        dvel = np.sqrt(rw_obs**2 + rn_obs**2 + (self.rr_obs - self.dist_sun)**2)

        # compute spacecraft velocity relative to each pixel
        self.v_obs = - (rw_obs * self.obs_vw + rn_obs * self.obs_vn + (self.rr_obs - self.dist_sun) * self.obs_vr) / dvel

        # cartesian coordinates
        x1 = rw_obs
        y1 = rn_obs * coscrlt + self.rr_obs * sincrlt
        z1 = -rn_obs * sincrlt + self.rr_obs * coscrlt

        hx = x1 * coscrln + z1 * sincrln
        hy = y1
        hz = -x1 * sincrln + z1 * coscrln

        # differential rotation
        omega = (np.pi/180.) * (1./86400.) * (14.713 - 2.396 * hy**2 - 1.787 * hy**4)

        # projection into rotating frame
        vx_rot = omega * hz * self.rsun_ref
        vy_rot = 0.0
        vz_rot = -omega * hx * self.rsun_ref

        v1 = coscrln * vx_rot - sincrln * vz_rot
        v2 = vy_rot
        v3 = sincrln * vx_rot + coscrln * vz_rot

        rot_vw = v1
        rot_vn = v2 * coscrlt - v3 * sincrlt
        rot_vr = v2 * sincrlt + v3 * coscrlt

        # compute differential rotation
        # TODO minus sign???? ask solaster people
        self.v_rot = (rw_obs * rot_vw + rn_obs * rot_vn + (self.rr_obs - self.dist_sun) * rot_vr) / dvel

        # correct the dopplergram by subtracting off velocities
        # self.image = self.image - v_rot - v_obs
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
            avg_int[i] = np.mean(ints[~np.isnan(ints)])

        # take averages in mu annuli to fit to
        mu_avgs = (mu_edge[1:] + mu_edge[0:-1]) / 2.0

        # set the initial guess parameters for optimization
        if self.is_continuum():
            p0 = [55925.8, 0.88, -0.23]
        elif self.is_filtergram():
            p0 = [1000, 0.081, 0.4998]
        else:
            return None

        # do the fit and divide out the LD profile
        popt, pcov = curve_fit(quad_darkening, mu_avgs, avg_int, p0=p0)
        self.ldark = quad_darkening(self.mu, *popt)
        self.iflat = self.image/self.ldark
        return None

    def rescale_to_hmi(self, hmi_image):
        assert self.is_filtergram()

        # rescale the image
        self.image, foot = reproject_interp((self.image, self.head), hmi_image.head)

        # TODO recalculate the mu???
        self.mu = hmi_image.mu

    def plot_image(self):
        if self.is_magnetogram():
            # get cmap
            cmap = plt.get_cmap("RdYlBu").copy()
            cmap.set_bad(color="black")

            # plot the sun
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            im = ax1.imshow(self.image, cmap=cmap, origin="lower", vmin=-4200, vmax=4200)
            cb = fig.colorbar(im)
            cb.set_label(r"${\rm Magnetic\ Field\ Strength\ (G)}$")
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm Corrected HMI\ LOS\ Magnetic\ Field}$")
            ax1.text(2650, 50, self.date_obs, fontsize=8, c="white")
            ax1.grid(False)
            fig.savefig("/Users/michael/Desktop/mag.pdf", bbox_inches="tight", dpi=500)
            plt.clf(); plt.close()
            return None

        elif self.is_dopplergram():
            # get cmap
            cmap = plt.get_cmap("seismic").copy()
            cmap.set_bad(color="black")

            # plot the sun
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            im = ax1.imshow(self.image - self.v_rot - self.v_obs, origin="lower", cmap=cmap, vmin=-2000, vmax=2000)
            cb = fig.colorbar(im)
            cb.set_label(r"${\rm LOS\ Velocity\ (km/s)}$")
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm Corrected\ HMI\ LOS\ Dopplergram}$")
            ax1.text(2650, 50, self.date_obs, fontsize=8, c="white")
            ax1.grid(False)
            fig.savefig("/Users/michael/Desktop/dop.pdf", bbox_inches="tight", dpi=500)
            plt.clf(); plt.close()
            return None

        elif self.is_continuum():
            # get cmap
            cmap = plt.get_cmap("afmhot").copy()
            cmap.set_bad(color="black")

            # plot the sun
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            im = ax1.imshow(self.iflat, cmap=cmap, origin="lower")#, vmin=20000)
            cb = fig.colorbar(im)
            cb.set_label(r"${\rm Relative\ Intensity}$")
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm Flattened\ HMI\ Continuum}$")
            ax1.text(2650, 50, self.date_obs, fontsize=8, c="white")
            ax1.grid(False)
            fig.savefig("/Users/michael/Desktop/con.pdf", bbox_inches="tight", dpi=500)
            plt.clf(); plt.close()
            return None

        elif self.is_filtergram():
            # get cmap
            cmap = plt.get_cmap("Purples_r").copy()
            cmap.set_bad(color="black")

            # plot the sun
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            im = ax1.imshow(self.image, cmap=cmap, origin="lower")#, vmin=20000)
            cb = fig.colorbar(im)
            cb.set_label(r"${\rm Relative\ Intensity}$")
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm Flattened\ AIA\ 1700\AA\ Filtergram}$")
            ax1.text(2650, 50, self.date_obs, fontsize=8, c="white")
            ax1.grid(False)
            fig.savefig("/Users/michael/Desktop/aia.pdf", bbox_inches="tight", dpi=500)
            plt.clf(); plt.close()
            return None

        else:
            return None

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

class SunMask:
    def __init__(self, con, mag, dop, aia):
        # check argument order/names are correct
        assert con.is_continuum()
        assert mag.is_magnetogram()
        assert dop.is_dopplergram()
        assert aia.is_filtergram()

        # calculate weights
        self.w_active, self.w_quiet = calculate_weights(mag)

        # identify regions
        self.identify_regions(con, mag, dop, aia)

        return None

    def identify_regions(self, con, mag, dop, aia):
        # calculate intensity thresholds for HMI and AIA
        con_thresh = 0.89 * np.nansum(con.iflat * self.w_quiet)/np.nansum(self.w_quiet)
        aia_thresh = np.nansum(aia.iflat * self.w_quiet)/np.nansum(self.w_quiet)

        # allocate memory for mask array
        self.regions = np.zeros(np.shape(con.image))

        # get thresholds for penumbrae, umbrae, and quiet sun
        ind1 = ((con.iflat <= con_thresh) & (con.iflat > (0.25 * con_thresh)))
        ind2 = (con.iflat <= (0.25 * con_thresh))
        ind3 = ((con.iflat > con_thresh))# & (aia.iflat < (1.3 * aia_thresh)))

        # set mask indices
        self.regions[ind1] = 1 # penumbrae
        self.regions[ind2] = 2 # umbrae
        self.regions[ind3] = 3 # quiet sun

        # bright region selection
        # ind4 = ((aia.iflat > (1.3 * aia_thresh)) & (mask != 1) & (mask != 2))
        # regions[ind4] = 4

        # make remaining regions quiet sun
        ind_rem = ((con.mu > 0.0) & (self.regions == 0))
        self.regions[ind_rem] = 3

        # set values beyond mu_thresh to nan
        self.regions[np.logical_or(con.mu <= con.mu_thresh, np.isnan(con.mu))] = np.nan

        return None

    def is_penumbra(self):
        return self.regions == 1

    def is_umbra(self):
        return self.regions == 2

    def is_quiet(self):
        return self.regions == 3

    def plot_image(self, date_obs):
        # get cmap
        cmap = colors.ListedColormap(["saddlebrown", "black", "orange"])
        cmap.set_bad(color="black")
        norm = colors.BoundaryNorm([0, 1, 2, 3], ncolors=cmap.N, clip=True)

        # plot the sun
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.imshow(self.regions - 0.5, cmap=cmap, norm=norm, origin="lower")
        cb = fig.colorbar(im, ticks=[0.5,1.5,2.5])
        cb.ax.set_yticklabels(["Penumbrae", "Umbrae", "Quiet Sun"])
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(r"${\rm SDO\ Identified\ Regions}$")
        ax1.text(2650, 50, date_obs, fontsize=8, c="white")
        ax1.grid(False)
        fig.savefig("/Users/michael/Desktop/mask.pdf", bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None
