# code copied/modified from https://github.com/samarth-kashyap/hmi-clean-ls
# original author: Samarth Ganesh Kashyap (g.samarth@tifr.res.in)
# original associated publication: https://arxiv.org/abs/2105.12055

import numpy as np
from math import pi
from pyshtools import legendre as pleg
from scipy.special import lpmv, lpmn, factorial, eval_legendre

def get_pleg_index(l, m):
    return int(l*(l+1)/2 + m)


def gen_leg(lmax, theta):
    cost = np.cos(theta)
    sint = np.sin(theta).reshape(1, theta.shape[0])

    maxIndex = int(lmax+1)
    ell = np.arange(maxIndex)
    norm = np.sqrt(ell*(ell+1)).reshape(maxIndex, 1)
    norm[norm == 0] = 1

    leg = np.zeros((maxIndex, theta.size))
    leg_d1 = np.zeros((maxIndex, theta.size))

    for i,z in enumerate(cost):
        leg[:, i], leg_d1[:, i] = pleg.PlBar_d1(lmax, z)
    return leg/np.sqrt(2)/norm, leg_d1 * (-sint)/np.sqrt(2)/norm

"""
Vectorized version of gen_leg()^
"""
def gen_leg_vec(lmax, theta):
    # build cos, sin and l arrays
    theta = theta.value * np.pi / 180.0
    cost = np.cos(theta)
    sint = np.sin(theta)
    ell = np.arange(lmax+1)

    # l(l+1)–scaling (avoid div 0 at l = 0)
    norm_l = np.sqrt(ell*(ell+1))
    norm_l[0] = 1

    # shtools “geodesy” normalization sqrt(2l+1)
    norm_sht = np.sqrt(2*ell + 1)

    # evaluate P_l(cost) for all l and all theta
    P = eval_legendre(ell[:,None], cost[None,:])

    # compute dP_l / dz via (1–z^2)P'_l(z) = l[P_(l-1)(z) – z P_l(z)] (recurrence relation)
    dP = np.zeros_like(P)
    z = cost[None,:]
    l = ell[1:,None]
    dP[1:,:] = l * (P[:-1,:] - z*P[1:,:]) / (1.0 - z*z)

    # apply sqrt(2l+1) normalization
    leg = P * norm_sht[:,None]
    leg_dz = dP * norm_sht[:,None]

    # final sqrt2 and sqrtl(l+1) scalings, and change of variable d theta = -sin theta dz
    leg_out =  leg / np.sqrt(2) / norm_l[:,None]
    leg_d1_out = -sint[None,:] * leg_dz / np.sqrt(2) / norm_l[:,None]
    return leg_out, leg_d1_out


def gen_leg_x(lmax, x):
    maxIndex = int(lmax+1)
    ell = np.arange(maxIndex)
    norm = np.sqrt(ell*(ell+1)).reshape(maxIndex, 1)
    norm[norm == 0] = 1

    leg = np.zeros((maxIndex, x.size))
    leg_d1 = np.zeros((maxIndex, x.size))

    for i,z in enumerate(x):
        leg[:, i], leg_d1[:, i] = pleg.PlBar_d1(lmax, z)
    return leg/np.sqrt(2)/norm, leg_d1/np.sqrt(2)/norm

"""
Vectorized version of gen_leg_x()^
"""
def gen_leg_x_vec(lmax, x):
    # switch to radians
    x = x.value * np.pi / 180.0

    # degree index l=0…lmax
    ell = np.arange(lmax+1)

    # l(l+1)–scaling (avoid div 0 at l = 0)
    norm_l = np.sqrt(ell*(ell+1))
    norm_l[0] = 1

    # shtools “geodesy” normalization sqrt(2l+1)
    norm_sht = np.sqrt(2*ell + 1)

    # evaluate all P_l(x)
    P = eval_legendre(ell[:,None], x[None,:])

    # compute dP_l/dx via (1−x^2) dP_l/dx = l [P_(l-1)(x) − x P_l(x)]
    dP = np.zeros_like(P)
    l = ell[1:,None]            # shape (lmax,1)
    z = x[None,:]
    dP[1:,:] = l * (P[:-1,:] - z*P[1:,:])/(1 - z*z)

    # apply sqrt(2l+1) normalization
    leg    = P    * norm_sht[:,None]
    leg_dx = dP   * norm_sht[:,None]

    # final sqrt2 and sqrtl(l+1) scalings
    leg_out    = leg    / np.sqrt(2) / norm_l[:,None]
    leg_d1_out = leg_dx / np.sqrt(2) / norm_l[:,None]
    return leg_out, leg_d1_out

def inv_SVD(A, svdlim):
    u, s, v = np.linalg.svd(A, full_matrices=False)
    sinv = s**-1
    sinv[sinv/sinv[0] > svdlim] = 0.0  # svdlim
    return np.dot(v.transpose().conjugate(),
                  np.dot(np.diag(sinv), u.transpose().conjugate()))
