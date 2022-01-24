import numpy as np
import sunpy as sp
from sunpy.net import Fido, attrs as a
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb
import glob
import warnings

def get_pix_grid(hdr):
    n_axis1  = hdr["NAXIS1"]                 # = 4096
    n_axis2  = hdr["NAXIS2"]                 # = 4096
    crpix1   = n_axis1 - hdr["CRPIX1"]
    crpix2   = n_axis2 - hdr["CRPIX2"]
    p_axis1 = np.arange(n_axis1) - crpix1
    p_axis2 = np.arange(n_axis2) - crpix2
    return np.meshgrid(p_axis1, p_axis2)

def get_geom_params(hdr):
    # get values from header
    cdelt1   = hdr["CDELT1"]
    cdelt2   = hdr["CDELT2"]
    dsun_obs = hdr["DSUN_OBS"]               # distance from observer to sun in m
    rsun_ref = hdr["RSUN_REF"]               # solar radius in m
    rsun_obs = hdr["RSUN_OBS"]               # apparent solar radius in arcsec

    # calculate solar radius in pixels on image
    dist_sun  = dsun_obs/rsun_ref               # distance to sun in solar radii
    foc       = 180. * 3600. / np.pi / cdelt1   # focal length in pixels
    p_sun     = foc * np.sqrt(1 - 1/dist_sun**2)/(dist_sun - 1/dist_sun)

    # make a mesh of distances
    px, py = get_pix_grid(hdr)
    dp = np.sqrt(px**2 + py**2)
    rr = (dist_sun * foc * dp - np.sqrt(foc**2 * dp**2 + dp**4 - dist_sun**2 * dp**4))/(foc**2 + dp**2)

    # make a dictionary
    soldict = {
        "dist_sun": dist_sun,
        "foc": foc,
        "p_sun": p_sun,
        "px": px,
        "py": py,
        "dp": dp,
        "rr": rr
    }
    return soldict

def get_mu_grid(soldict):
    # now do the trig
    cos_alpha = np.sqrt(1.0 - soldict["rr"]**2)
    sin_alpha = soldict["rr"]
    cos_theta = (soldict["dist_sun"] - cos_alpha) / np.sqrt(soldict["rr"]**2 + (soldict["dist_sun"] - cos_alpha)**2)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    return cos_alpha * cos_theta - sin_alpha * sin_theta

def get_differential_rotation(hdr, soldict):
    # get values from header
    crota2   = hdr["CROTA2"] * np.pi/180.0
    crln_obs = hdr["CRLN_OBS"] * np.pi/180.0 # Carrington Longitude of observer
    crlt_obs = hdr["CRLT_OBS"] * np.pi/180.0 # Heliographic latitude of observer
    rsun_ref  = hdr["RSUN_REF"]
    rsun_obs  = hdr["RSUN_OBS"]

    # transforms
    coscrlt  = np.cos(crlt_obs)
    sincrlt  = np.sin(crlt_obs)
    coscrln  = np.cos(crln_obs)
    sincrln  = np.sin(crln_obs)
    sincrota = np.sin(crota2)
    coscrota = np.cos(crota2)

    # do some math
    dw = soldict["py"] * sincrota + soldict["px"] * coscrota
    dn = soldict["py"] * coscrota - soldict["px"] * sincrota
    rw_obs = soldict["rr"] * dw/soldict["dp"]
    rn_obs = soldict["rr"] * dn/soldict["dp"]
    rr_obs = np.sqrt(1.0 - soldict["rr"]**2)

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
    vx_rot = omega * hz * rsun_ref
    vy_rot = 0
    vz_rot = -omega * hx * rsun_ref
    v1 = coscrln * vx_rot - sincrln * vz_rot
    v2 = vy_rot
    v3 = sincrln * vx_rot + coscrln * vz_rot
    rot_vw = v1
    rot_vn = v2 * coscrlt - v3 * sincrlt
    rot_vr = v2 * sincrlt + v3 * coscrlt

    # calculate grid
    num = rw_obs * rot_vw + rn_obs * rot_vn + (rr_obs - soldict["dist_sun"]) * rot_vr
    den = np.sqrt(rw_obs**2 + rn_obs**2 + (rr_obs - soldict["dist_sun"])**2)
    return num/den

def get_observer_velocity(hdr, soldict):
    # get values from header
    crota2   = hdr["CROTA2"] * np.pi/180.0
    obs_vr   = hdr["OBS_VR"]                 # radial component of spacecraft velocity
    obs_vw   = hdr["OBS_VW"]                 # westward  ''
    obs_vn   = hdr["OBS_VN"]                 # northward ''

    # transforms
    sincrota = np.sin(crota2)
    coscrota = np.cos(crota2)

    # do some math
    dw = soldict["py"] * sincrota + soldict["px"] * coscrota
    dn = soldict["py"] * coscrota - soldict["px"] * sincrota
    rw_obs = soldict["rr"] * dw/soldict["dp"]
    rn_obs = soldict["rr"] * dn/soldict["dp"]
    rr_obs = np.sqrt(1.0 - soldict["rr"]**2)

    # get observer velocity
    num = -(rw_obs * obs_vw + rn_obs * obs_vn + (rr_obs - soldict["dist_sun"]) * obs_vr)
    den = np.sqrt(rw_obs**2 + rn_obs**2 + (rr_obs - soldict["dist_sun"])**2)
    return num/den
