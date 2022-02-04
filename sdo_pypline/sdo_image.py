import numpy as np
from .sdo_io import *
from .limbdark import *

import pdb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from reproject import reproject_interp
from astropy.io.fits import Header, FITSFixedWarning
warnings.filterwarnings("ignore", category=FITSFixedWarning)

class HMI_Image:
    def __init__(self, file):
        # read the image file
        self.file = file
        head = Header(read_header(file))

        # parse out the header
        self.naxis1 = head["NAXIS1"]
        self.naxis2 = head["NAXIS2"]
        self.wavelength = head["WAVELNTH"]
        self.crpix1 = head["CRPIX1"]
        self.crpix2 = head["CRPIX2"]
        self.crval1 = head["CRVAL1"]
        self.crval2 = head["CRVAL2"]
        self.cdelt1 = head["CDELT1"]
        self.cdelt2 = head["CDELT2"]
        self.crota2 = head["CROTA2"] * np.pi/180.0
        self.dsun_obs = head["DSUN_OBS"]
        self.dsun_ref = head["DSUN_REF"]
        self.rsun_obs = head["RSUN_OBS"]
        self.rsun_ref = head["RSUN_REF"]
        self.crln_obs = head["CRLN_OBS"] * np.pi/180.0
        self.crlt_obs = head["CRLT_OBS"] * np.pi/180.0
        self.car_rot = head["CAR_ROT"]
        self.obs_vr = head["OBS_VR"]
        self.obs_vw = head["OBS_VW"]
        self.obs_vn = head["OBS_VN"]
        self.rsun_obs = head["RSUN_OBS"]
        self.date_obs = head["DATE-OBS"]
        self.content = head["CONTENT"]

        # now set head attribute
        self.head = head

        # get mesh of distances and pixels
        self.dist_sun = self.dsun_obs/self.rsun_ref
        self.focal_len = 180. * 3600. / np.pi / self.cdelt1
        self.p_sun = self.focal_len * np.sqrt(1. - 1./self.dist_sun**2.0)/ \
                     (self.dist_sun - 1./self.dist_sun)

        # make a mesh of distances
        paxis1 = np.arange(self.naxis1) - self.crpix1
        paxis2 = np.arange(self.naxis2) - self.crpix2
        self.px, self.py = np.meshgrid(paxis1, paxis2)
        self.dp = np.sqrt(self.px**2.0 + self.py**2.0)
        self.rr = (self.dist_sun * self.focal_len * self.dp - \
                   np.sqrt(self.focal_len**2 * self.dp**2 + self.dp**4 - self.dist_sun**2 * self.dp**4))/ \
                   (self.focal_len**2 + self.dp**2)

        # calculate grid of mus
        cos_alpha = np.sqrt(1.0 - self.rr**2)
        sin_alpha = self.rr
        cos_theta = (self.dist_sun - cos_alpha) / np.sqrt(self.rr**2 + (self.dist_sun - cos_alpha)**2)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        self.mu = cos_alpha * cos_theta - sin_alpha * sin_theta

        # read in the data
        self.image = read_data(file)

    def is_magnetogram(self):
        return self.content == "MAGNETOGRAM"

    def is_dopplergram(self):
        return self.content == "DOPPLERGRAM"

    def is_continuum(self):
        return self.content == "CONTINUUM INTENSITY"

    def mask_low_mu(self, mu_thresh):
        self.image[np.logical_or(self.mu <= mu_thresh, np.isnan(self.mu))] = np.nan

    def calc_differential_rot(self):
        # geometric quantities
        coscrlt  = np.cos(self.crlt_obs)
        sincrlt  = np.sin(self.crlt_obs)
        coscrln  = np.cos(self.crln_obs)
        sincrln  = np.sin(self.crln_obs)
        sincrota = np.sin(self.crota2)
        coscrota = np.cos(self.crota2)

        # do the trig
        dw = self.py * sincrota + self.px * coscrota
        dn = self.py * coscrota - self.px * sincrota
        rw_obs = self.rr * dw/self.dp
        rn_obs = self.rr * dn/self.dp
        rr_obs = np.sqrt(1.0 - self.rr**2)

        # geometry
        x1 = rw_obs
        y1 = rn_obs * coscrlt + rr_obs * sincrlt
        z1 = -rn_obs * sincrlt + rr_obs * coscrlt
        hx = x1 * coscrln + z1 * sincrln
        hy = y1
        hz = -x1 * sincrln + z1 * coscrln

        # differential rotation
        omega  = (1./86400.) * (np.pi/180.) * (14.713 - 2.396 * hy**2 - 1.787 * hy**4)

        # projection
        vx_rot = omega * hz * self.rsun_ref
        vy_rot = 0
        vz_rot = -omega * hx * self.rsun_ref
        v1 = coscrln * vx_rot - sincrln * vz_rot
        v2 = vy_rot
        v3 = sincrln * vx_rot + coscrln * vz_rot
        rot_vw = v1
        rot_vn = v2 * coscrlt - v3 * sincrlt
        rot_vr = v2 * sincrlt + v3 * coscrlt

        # calculate grid
        num = rw_obs * rot_vw + rn_obs * rot_vn + (rr_obs - self.dist_sun) * rot_vr
        den = np.sqrt(rw_obs**2 + rn_obs**2 + (rr_obs - self.dist_sun)**2)
        return num/den

    def calc_observer_vel(self):
        # transforms
        sincrota = np.sin(self.crota2)
        coscrota = np.cos(self.crota2)

        # do some math
        dw = self.py * sincrota + self.px * coscrota
        dn = self.py * coscrota - self.px * sincrota
        rw_obs = self.rr * dw/self.dp
        rn_obs = self.rr * dn/self.dp
        rr_obs = np.sqrt(1.0 - self.rr**2)

        # get observer velocity
        num = -(rw_obs * self.obs_vw + rn_obs * self.obs_vn + (rr_obs - self.dist_sun) * self.obs_vr)
        den = np.sqrt(rw_obs**2 + rn_obs**2 + (rr_obs - self.dist_sun)**2)
        return num/den

    def correct_magnetogram(self):
        assert self.is_magnetogram()
        self.image /= self.mu

    def correct_dopplergram(self):
        # only do the correction if its a dopplergram
        assert self.is_dopplergram()

        # compute the velocities
        # TOOD: check on the rotation???
        # obs_vel = self.calc_observer_vel()
        # sun_vel = self.calc_differential_rot()
        obs_vel = np.rot90(np.flip(self.calc_observer_vel(), axis=1), 2)
        sun_vel = np.rot90(np.flip(self.calc_differential_rot(), axis=1), 2)

        # update image data
        self.image = self.image - obs_vel - sun_vel
        return self.image

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


class AIA_Image:
    def __init__(self, file):
        # read the image file
        self.file = file
        head = Header(read_header(file))

        # parse out the header
        self.naxis1 = head["NAXIS1"]
        self.naxis2 = head["NAXIS2"]
        self.wavelength = head["WAVELNTH"]
        self.crpix1 = head["CRPIX1"]
        self.crpix2 = head["CRPIX2"]
        self.crval1 = head["CRVAL1"]
        self.crval2 = head["CRVAL2"]
        self.cdelt1 = head["CDELT1"]
        self.cdelt2 = head["CDELT2"]
        self.crota2 = head["CROTA2"] * np.pi/180.0
        self.dsun_obs = head["DSUN_OBS"]
        self.dsun_ref = head["DSUN_REF"]
        self.rsun_obs = head["RSUN_OBS"]
        self.rsun_ref = head["RSUN_REF"]
        self.crln_obs = head["CRLN_OBS"] * np.pi/180.0
        self.crlt_obs = head["CRLT_OBS"] * np.pi/180.0
        self.car_rot = head["CAR_ROT"]
        self.obs_vr = head["OBS_VR"]
        self.obs_vw = head["OBS_VW"]
        self.obs_vn = head["OBS_VN"]
        self.rsun_obs = head["RSUN_OBS"]
        self.date_obs = head["DATE-OBS"]
        self.content = "FILTERGRAM"

        # now set head attribute
        self.head = head

        # get mesh of distances and pixels
        self.dist_sun = self.dsun_obs/self.rsun_ref
        self.focal_len = 180. * 3600. / np.pi / self.cdelt1
        self.p_sun = self.focal_len * np.sqrt(1. - 1./self.dist_sun**2.0)/ \
                     (self.dist_sun - 1./self.dist_sun)

        # make a mesh of distances
        paxis1 = np.arange(self.naxis1) - self.crpix1
        paxis2 = np.arange(self.naxis2) - self.crpix2
        self.px, self.py = np.meshgrid(paxis1, paxis2)
        self.dp = np.sqrt(self.px**2.0 + self.py**2.0)
        self.rr = (self.dist_sun * self.focal_len * self.dp - \
                   np.sqrt(self.focal_len**2 * self.dp**2 + self.dp**4 - self.dist_sun**2 * self.dp**4))/ \
                   (self.focal_len**2 + self.dp**2)

        # calculate grid of mus
        cos_alpha = np.sqrt(1.0 - self.rr**2)
        sin_alpha = self.rr
        cos_theta = (self.dist_sun - cos_alpha) / np.sqrt(self.rr**2 + (self.dist_sun - cos_alpha)**2)
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        self.mu = cos_alpha * cos_theta - sin_alpha * sin_theta

        # read in the data
        self.image = read_data(file).astype(float)

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
        self.image, foot = reproject_interp((self.image, self.head), hmi_image.head)
