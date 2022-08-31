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
import argparse, warnings
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
    parser.add_argument('end', type=str, help='ending date formatted as YYYY/MM/DD')
    parser.add_argument('sample', type=int, help='cadence of sampling in hours')

    # parse the command line arguments
    args = parser.parse_args()
    outdir = args.outdir
    start = args.start
    end = args.end
    sample = args.sample

    # set time attributes for search
    start += 'T00:00:00'
    end += 'T24:00:00'
    trange = a.Time(start, end)
    sample = a.Sample(sample * u.hour)
    provider = a.Provider("JSOC")

    # set attributes for HMI query
    instr1 = a.Instrument.hmi
    physobs = (a.Physobs.los_velocity | a.Physobs.los_magnetic_field | a.Physobs.intensity)

    # set attributes for AIA query
    level = a.Level(1)
    instr2 = a.Instrument.aia
    wavelength = a.Wavelength(1700. * u.AA)

    # get query for HMI and download data
    qr1 = Fido.search(trange, instr1, physobs, provider, sample)
    hmi_files = Fido.fetch(qr1, path=outdir, overwrite=False, progress=True)

    # get query for AIA and download data
    qr2 = Fido.search(trange, instr2, wavelength, level, provider, sample)
    aia_files = Fido.fetch(qr2, path=outdir, overwrite=False, progress=True)

if __name__ == '__main__':
    main()
