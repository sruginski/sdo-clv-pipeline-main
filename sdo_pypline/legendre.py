# code copied/modified from https://github.com/samarth-kashyap/hmi-clean-ls
# original author: Samarth Ganesh Kashyap (g.samarth@tifr.res.in)
# original associated publication: https://arxiv.org/abs/2105.12055

import numpy as np
from math import pi
from pyshtools import legendre as pleg

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

    count = 0
    for z in cost:
        leg[:, count], leg_d1[:, count] = pleg.PlBar_d1(lmax, z)
        count += 1
    return leg/np.sqrt(2)/norm, leg_d1 * (-sint)/np.sqrt(2)/norm


def gen_leg_x(lmax, x):
    maxIndex = int(lmax+1)
    ell = np.arange(maxIndex)
    norm = np.sqrt(ell*(ell+1)).reshape(maxIndex, 1)
    norm[norm == 0] = 1

    leg = np.zeros((maxIndex, x.size))
    leg_d1 = np.zeros((maxIndex, x.size))

    count = 0
    for z in x:
        leg[:, count], leg_d1[:, count] = pleg.PlBar_d1(lmax, z)
        count += 1
    return leg/np.sqrt(2)/norm, leg_d1/np.sqrt(2)/norm


def inv_SVD(A, svdlim):
    u, s, v = np.linalg.svd(A, full_matrices=False)
    sinv = s**-1
    sinv[sinv/sinv[0] > svdlim] = 0.0  # svdlim
    return np.dot(v.transpose().conjugate(),
                  np.dot(np.diag(sinv), u.transpose().conjugate()))
