# code copied/modified from https://github.com/samarth-kashyap/hmi-clean-ls
# original author: Samarth Ganesh Kashyap (g.samarth@tifr.res.in)
# original associated publication: https://arxiv.org/abs/2105.12055

import numpy as np
from math import pi
from pyshtools import legendre as pleg
from scipy.special import lpmv

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
    # array of degrees 0…lmax
    x = np.cos(theta)
    L = np.arange(lmax + 1)

    # P[k, i] = P_k(x[i])
    P = np.vstack([lpmv(0, k, x) for k in L])

    # compute dP/dx via three-term recurrence
    dP = np.zeros_like(P)
    for k in range(1, lmax + 1):
        dP[k] = (k * x * P[k] - k * P[k - 1]) / (x*x - 1)

    # convert to dP/dθ = -sin(theta) * dP/dx
    dP_dtheta = -np.sin(theta)[None, :] * dP

    # normalization factor sqrt(l*(l+1)), with entry 0 set to 1
    norm = np.sqrt(L*(L + 1))[:, None]
    norm[0] = 1

    # apply 1/sqrt(2) and norm
    Pn = P / np.sqrt(2) / norm
    dPn = dP_dtheta / np.sqrt(2) / norm
    return Pn, dPn

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
    # array of degrees 0…lmax
    L = np.arange(lmax + 1)

    # P[k,i] = P_k(x[i])
    P = np.vstack([lpmv(0, k, x) for k in L])

    # dP/dx via recurrence
    dP = np.zeros_like(P)
    for k in range(1, lmax + 1):
        dP[k] = (k * x * P[k] - k * P[k-1]) / (x*x - 1)

    # normalization factor sqrt(l*(l+1)), with entry 0 set to 1
    norm = np.sqrt(L*(L + 1))[:, None]
    norm[0] = 1

    # apply 1/sqrt(2) and norm
    Pn  = P / np.sqrt(2) / norm
    dPn = dP / np.sqrt(2) / norm
    return Pn, dPn


def inv_SVD(A, svdlim):
    u, s, v = np.linalg.svd(A, full_matrices=False)
    sinv = s**-1
    sinv[sinv/sinv[0] > svdlim] = 0.0  # svdlim
    return np.dot(v.transpose().conjugate(),
                  np.dot(np.diag(sinv), u.transpose().conjugate()))
