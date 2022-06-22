# necessary modules
import numpy as np
import sunpy as sp
import datetime as dt
import re, pdb, csv, glob
from astropy.io import fits
from astropy.time import Time
from os.path import exists, split, isdir, getsize

# read headers and data
def read_header(file):
    # return fits.getheader(file, 1, output_verify="silentfix")
    with fits.open(file) as hdu_list:
        hdu_list.verify("silentfix")
        header = hdu_list[1].header
    return header

def read_data(file):
    with fits.open(file) as hdu_list:
        hdu_list.verify("silentfix")
        data = hdu_list[1].data.astype(float)
    return data

# function to glob the input data
def find_data(indir):
    # find the data
    con_files, con_dates = sort_data(glob.glob(indir + "*hmi.Ic*.fits"))
    mag_files, mag_dates = sort_data(glob.glob(indir + "*hmi.M*.fits"))
    dop_files, dop_dates = sort_data(glob.glob(indir + "*hmi.v*.fits"))
    aia_files, aia_dates = sort_data(glob.glob(indir + "*aia*.fits"))

    # find datetimes that are in *all* lists
    common_dates = list(set.intersection(*map(set, [con_dates, mag_dates, dop_dates, aia_dates])))

    # remove epochs that are missing in any data set from all data sets
    con_files = [con_files[idx] for idx, date in enumerate(con_dates) if date in common_dates]
    mag_files = [mag_files[idx] for idx, date in enumerate(mag_dates) if date in common_dates]
    dop_files = [dop_files[idx] for idx, date in enumerate(dop_dates) if date in common_dates]
    aia_files = [aia_files[idx] for idx, date in enumerate(aia_dates) if date in common_dates]
    return con_files, mag_files, dop_files, aia_files

def sort_data(f_list):
    # get the dates and indices to sort
    dates = get_dates(f_list)
    inds = np.argsort(dates)
    return [f_list[i] for i in inds], [dates[i] for i in inds]

def get_date(f):
    if "aia" in f:
        s = re.search(r'\d{4}-\d{2}-\d{2}T\d{2}\d{2}\d{2}', f)
    else:
        s = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}', f)
        if s is None:
            s = re.search(r'\d{4}.\d{2}.\d{2}_\d{2}_\d{2}_\d{2}', f)


    # standardize time formats
    s = s.group()
    s = s.replace("t", "_")
    s = s.replace("T", "_")
    s = s.replace("-", "")
    return round_time(date=dt.datetime.strptime(s, "%Y%m%d_%H%M%S"))

def get_dates(files):
    date = []
    for f in files:
        # make it a datetime
        date.append(get_date(f))

    return date

def round_time(date=None, round_to=3600):
   """Round a datetime object to any time lapse in seconds
   dt : datetime.datetime object, default now.
   round_to : Closest number of seconds to round to, default 1 hour.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   if date == None : date = datetime.datetime.now()
   seconds = (date.replace(tzinfo=None) - date.min).seconds
   rounding = (seconds+round_to/2) // round_to * round_to
   return date + dt.timedelta(0,rounding-seconds,-date.microsecond)

def truncate_output_file(fname):
    # truncate the file if it does exist
    if exists(fname):
        with open(fname, "w") as f:
            f.truncate()
    return None

def create_output_file(fname):
    if not exists(fname):
        with open(fname, "w") as f:
            pass

def find_last_date(fname):
    with open(fname, "r") as f:
        for line in f:
            pass
        mjd_str = line.split(",")[0]
    return mjd_str

def remove_line_by_mjd(mjd_str, fname):
    # find the lines that don't include this mjd
    lines = []
    with open(fname, "r") as f:
        reader = csv.reader(f)
        for idx, row in enumerate(reader):
            if (idx == 0) | (mjd_str not in row):
                lines.append(row)

    # wipe the file
    truncate_file(fname)

    # write good lines to file
    with open(fname, "w") as f:
        writer = csv.writer(f)
        writer.writerows(lines)

    return None

def create_file(fname, header):
    with open(fname, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
    return None

def write_vels(fname, mjd, ffactor, Bobs, pen_frac, umb_frac, quiet_frac, plage_frac, vels):
    # create the file if it doesn't exist or if it's empty
    if (((not exists(fname)) or getsize(fname) == 0) and isdir(split(fname)[0])):
        create_file(fname, ["mjd", "ffactor", "Bobs", "pen_frac", "umb_frac", \
                            "quiet_frac", "plage_frac", "v_hat", \
                            "v_phot", "v_quiet", "v_conv"])

    # write the vels
    with open(fname, "a") as f:
        writer = csv.writer(f)
        writer.writerow(np.concatenate(([mjd, ffactor, Bobs, pen_frac, umb_frac, quiet_frac, plage_frac], [v for v in vels])))

    return None

def write_vels_by_region(fname, mjd, region, lo_mu, hi_mu, vels):
    # create the file if it doesn't exist or if it's empty
    if (((not exists(fname)) or getsize(fname) == 0) and isdir(split(fname)[0])):
         create_file(fname, ["mjd", "region", "lo_mu", "hi_mu", \
                             "v_hat", "v_phot", "v_quiet", "v_conv"])

    # write the vels
    with open(fname, "a") as f:
        writer = csv.writer(f)
        writer.writerow(np.concatenate(([mjd], [region], [lo_mu], [hi_mu], [v for v in vels])))

    return None
