import numpy as np
from .sdo_io import *
from .limbdark import *

import pdb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from reproject import reproject_interp
from astropy.io.fits import Header#, FITSFixedWarning
# warnings.filterwarnings("ignore", category=FITSFixedWarning)


class HMI_Image:
    def __init__(self, file):
        # set attributes from the heaader
        self.parse_header(file)

        # read in the data
        self.image = read_data(file)

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
        self.rr_obs = np.sqrt(1.0 - self.rr**2)

        # calculate grid of mus
        cos_alpha = self.rr_obs
        sin_alpha = self.rr
        cos_theta = (self.dist_sun - cos_alpha) / np.sqrt(self.rr**2 + (self.dist_sun - cos_alpha)**2)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        self.mu = np.real(cos_alpha * cos_theta - sin_alpha * sin_theta)

    def parse_header(self, file):
        # read the header
        head = Header(read_header(file))
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
        self.content = head["CONTENT"]

    def is_magnetogram(self):
        return self.content == "MAGNETOGRAM"

    def is_dopplergram(self):
        return self.content == "DOPPLERGRAM"

    def is_continuum(self):
        return self.content == "CONTINUUM INTENSITY"

    def mask_low_mu(self, mu_thresh):
        self.image[np.logical_or(self.mu <= mu_thresh, np.isnan(self.mu))] = np.nan
        return None

    def correct_magnetogram(self):
        assert self.is_magnetogram()
        self.image /= self.mu
        return None

    def correct_dopplergram(self):
        # only do the correction if its a dopplergram
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
        v_obs = - (rw_obs * self.obs_vw + rn_obs * self.obs_vn + (self.rr_obs - self.dist_sun) * self.obs_vr) / dvel

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
        v_rot = (rw_obs * rot_vw + rn_obs * rot_vn + (self.rr_obs - self.dist_sun) * rot_vr) / dvel

        # correct the dopplergram by subtracting off velocities
        self.image = self.image - v_rot - v_obs
        return None

    def correct_limb_darkening(self, mu_lim=0.1, num_mu=50):
        # only do the correction if its a dopplergram
        assert self.is_continuum()

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

        # fit the data
        mu_avgs = (mu_edge[1:] + mu_edge[0:-1]) / 2.0
        p0  = [55925.8, 0.88, -0.23]
        popt, pcov = curve_fit(quad_darkening, mu_avgs, avg_int, p0=p0)
        self.image /= quad_darkening(self.mu, *popt)
        return None

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
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm HMI\ LOS\ Magnetic\ Field\ Strength}$")
            ax1.text(2700, 50, self.date_obs, fontsize=8, c="white")
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
            im = ax1.imshow(self.image, origin="lower", cmap=cmap, vmin=-2000, vmax=2000)
            cb = fig.colorbar(im)
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm HMI\ LOS\ Doppler\ Velocity}$")
            ax1.text(2700, 50, self.date_obs, fontsize=8, c="white")
            ax1.grid(False)
            fig.savefig("/Users/michael/Desktop/dop.pdf", bbox_inches="tight", dpi=500)
            plt.clf(); plt.close()
            return None

        elif self.is_continuum():
            # get cmap
            cmap = plt.get_cmap("afmhot").copy()
            cmap.set_bad(color="white")

            # plot the sun
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            im = ax1.imshow(self.image, cmap=cmap, origin="lower")#, vmin=20000)
            cb = fig.colorbar(im)
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.set_title(r"${\rm HMI\ Continuum\ Intensity}$")
            ax1.text(2700, 50, self.date_obs, fontsize=8)
            ax1.grid(False)
            fig.savefig("/Users/michael/Desktop/con.pdf", bbox_inches="tight", dpi=500)
            plt.clf(); plt.close()
            return None

        else:
            return None


class AIA_Image:
    def __init__(self, file):
        # set attributes from the heaader
        self.parse_header(file)

        # read in the data
        self.image = read_data(file)

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
        self.rr_obs = np.sqrt(1.0 - self.rr**2)

        # calculate grid of mus
        cos_alpha = self.rr_obs
        sin_alpha = self.rr
        cos_theta = (self.dist_sun - cos_alpha) / np.sqrt(self.rr**2 + (self.dist_sun - cos_alpha)**2)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        self.mu = np.real(cos_alpha * cos_theta - sin_alpha * sin_theta)

    def parse_header(self, file):
        # read the header
        head = Header(read_header(file))
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
        self.content = "FILTERGRAM"


    def correct_limb_darkening(self, mu_lim=0.1, num_mu=50):
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

        # fit the data
        mu_avgs = (mu_edge[1:] + mu_edge[0:-1]) / 2.0
        p0  = [1000, 0.081, 0.4998]
        popt, pcov = curve_fit(quad_darkening, mu_avgs, avg_int, p0=p0)
        self.image /= quad_darkening(self.mu, *popt)

    def rescale_to_hmi(self, hmi_image):
        # rescale the image
        self.image, foot = reproject_interp((self.image, self.head), hmi_image.head)
        self.mu = hmi_image.mu

    def plot_image(self):
        # get cmap
        cmap = plt.get_cmap("Purples_r").copy()
        cmap.set_bad(color="black")

        # plot the sun
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im = ax1.imshow(self.image, cmap=cmap, origin="lower")#, vmin=20000)
        cb = fig.colorbar(im)
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(r"${\rm AIA\ Filtergram}$")
        ax1.text(2700, 50, self.date_obs, fontsize=8, c="white")
        ax1.grid(False)
        fig.savefig("/Users/michael/Desktop/aia.pdf", bbox_inches="tight", dpi=500)
        plt.clf(); plt.close()
        return None
