# necessary modules
import numpy as np
import sunpy as sp
import datetime as dt
import re, pdb, csv, glob
from astropy.io import fits
from os.path import exists, split, isdir, getsize

# read headers and data
def read_header(file):
    return fits.getheader(file, 1, output_verify="silentfix")

def read_data(file):
    return fits.getdata(file, 1, output_verify="silentfix").astype(float)

# function to glob the input data
def find_data(indir):
    # find the data
    con_files = sort_data(glob.glob(indir + "*hmi_ic*.fits"))
    mag_files = sort_data(glob.glob(indir + "*hmi_m*.fits"))
    dop_files = sort_data(glob.glob(indir + "*hmi_v*.fits"))
    aia_files = sort_data(glob.glob(indir + "*aia*.fits"))
    return con_files, mag_files, dop_files, aia_files

def get_date(f):
    s = re.search(r'\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}', f)
    if s is None:
        s = re.search(r'\d{4}_\d{2}_\d{2}t\d{2}_\d{2}_\d{2}', f)

    # replace any t's with underscore
    s = s.group()
    s = s.replace("t", "_")
    return dt.datetime.strptime(s, "%Y_%m_%d_%H_%M_%S")

def get_dates(files):
    date = []
    for f in files:
        # make it a datetime
        date.append(get_date(f))

    return date

def sort_data(f_list):
    # get the dates and indices to sort
    dates = get_dates(f_list)
    inds = np.argsort(dates)
    return [f_list[i] for i in inds]

def truncate_file(fname):
    # truncate the file if it does exist
    if exists(fname):
        with open(fname, "w") as f:
            f.truncate()
    return fname

def create_file(fname, header):
    with open(fname, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
    return None

def write_vels(fname, mjd, ffactor, pen_frac, umb_frac, quiet_frac, plage_frac, vels):
    # create the file if it doesn't exist or if it's empty
    if (((not exists(fname)) or getsize(fname) == 0) and isdir(split(fname)[0])):
        create_file(fname, ["mjd", "ffactor", "pen_frac", "umb_frac", \
                            "quiet_frac", "plage_frac", "v_hat", \
                            "v_phot", "v_quiet", "v_conv"])

    # write the vels
    with open(fname, "a") as f:
        writer = csv.writer(f)
        writer.writerow(np.concatenate(([mjd, ffactor, pen_frac, umb_frac, quiet_frac, plage_frac], [v for v in vels])))

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
