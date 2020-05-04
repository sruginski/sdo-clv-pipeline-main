#==============================================================================
# Author: Michael Palumbo
# Date: January 2019
# Purpose: Download SDO data from JSOC
#==============================================================================
import numpy as np
import pandas as pd
import datetime as dt
import astropy.units as u
import os, sys, pdb, time
import argparse
import warnings
from sunpy.net import Fido, attrs as a
from sunpy.time import TimeRange
from astropy.units.quantity import AstropyDeprecationWarning

def main():
    # supress warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=AstropyDeprecationWarning)

    # initialize argparser
    parser = argparse.ArgumentParser(description='Download SDO data from JSOC',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('outdir', type=str, help='full directory path for file writeout')
    parser.add_argument('start', type=str, help='starting date formatted as YYYY/MM/DD')
    parser.add_argument('nstep', type=int, help='number of four hour steps from start time')

    # parse the command line arguments
    args = parser.parse_args()
    outdir = args.outdir
    start = args.start
    nstep = args.nstep

    # get start time
    start += 'T00:00:00'
    delta = dt.datetime.strptime(start, '%Y/%m/%dT%H:%M:%S') + dt.timedelta(hours=nstep * 4)

    # set search range for time
    range1 = TimeRange(delta, 120 * u.second)
    range2 = TimeRange(delta, 20 * u.second)

    # create queries for data (HMI then AIA)
    qr1 = Fido.search(a.Time(range1.start, range1.end),
                      a.jsoc.Notify('mlp95@psu.edu'),
                      a.jsoc.Series('hmi.M_720s') |
                      a.jsoc.Series('hmi.V_720s') |
                      a.jsoc.Series('hmi.Ic_720s'))
    qr2 = Fido.search(a.Time(range2.start, range2.end),
                      a.jsoc.Series('aia.lev1_uv_24s'),
                      a.jsoc.Wavelength(1700 * u.AA),
                      a.jsoc.Notify('mlp95@psu.edu'))

    # download the data from the queries
    req1 = Fido.fetch(qr1, path=outdir, overwrite=True, progress=True)
    req2 = Fido.fetch(qr2, path=outdir, overwrite=True, progress=True)
    time.sleep(6)

if __name__ == '__main__':
    main()
