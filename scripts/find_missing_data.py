import numpy as np
import datetime as dt
import astropy.units as u
import matplotlib.pyplot as plt
import re, pdb, csv, glob, argparse
from astropy.time import Time
from os.path import exists, split, isdir
from sunpy.net import Fido, attrs as a
from sunpy.time import TimeRange
from astropy.units.quantity import AstropyDeprecationWarning


from sdo_pypline.sdo_io import *

def main():
    # data download params
    sample1 = 6
    sample2 = a.Sample(6 * u.hour)
    provider = a.Provider("JSOC")
    level = a.Level(1)
    instr1 = a.Instrument.hmi
    instr2 = a.Instrument.aia
    physobs = (a.Physobs.los_velocity | a.Physobs.los_magnetic_field | a.Physobs.intensity)
    wavelength = a.Wavelength(1700. * u.AA)

    # find the data
    # TODO file names
    indir = "/Users/michael/Desktop/sdo_data/"
    con_files, con_dates = sort_data(glob.glob(indir + "*hmi_ic*.fits"))
    mag_files, mag_dates = sort_data(glob.glob(indir + "*hmi_m*.fits"))
    dop_files, dop_dates = sort_data(glob.glob(indir + "*hmi_v*.fits"))
    aia_files, aia_dates = sort_data(glob.glob(indir + "*aia*.fits"))

    # get earliest and latest date
    min_date = min(map(min, (con_dates, mag_dates, dop_dates, aia_dates)))
    max_date = max(map(max, (con_dates, mag_dates, dop_dates, aia_dates)))

    # get array of all datetimes that should be there
    dts = np.arange(min_date, max_date, dt.timedelta(hours=sample1)).astype(dt.datetime)

    for t in dts:
        if any([(t not in con_dates), (t not in mag_dates), (t not in dop_dates), (t not in aia_dates)]):
            # get t_range
            trange = a.Time(t, t + dt.timedelta(hours=1))

            # get query for HMI and download data
            qr1 = Fido.search(trange, instr1, physobs, provider, sample2)
            hmi_files = Fido.fetch(qr1, path=indir, overwrite=False, progress=False)

            # get query for AIA and download data
            qr2 = Fido.search(trange, instr2, wavelength, level, provider, sample2)
            aia_files = Fido.fetch(qr2, path=indir, overwrite=False, progress=False)
        else:
            continue

    pdb.set_trace()

    return None


if __name__ == "__main__":
    main()
