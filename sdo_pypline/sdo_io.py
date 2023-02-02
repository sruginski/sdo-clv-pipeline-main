# necessary modules
import numpy as np
import sunpy as sp
import datetime as dt
import os, re, pdb, csv, glob
from astropy.io import fits
from astropy.time import Time
from os.path import exists, split, isdir, getsize, splitext

from .paths import root

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
def find_data(indir, globexp=""):
    # find the data
    con_files, con_dates = sort_data(glob.glob(indir + "*hmi*" + globexp + "*con*.fits"))
    mag_files, mag_dates = sort_data(glob.glob(indir + "*hmi*" + globexp + "*mag*.fits"))
    dop_files, dop_dates = sort_data(glob.glob(indir + "*hmi*" + globexp + "*dop*.fits"))
    aia_files, aia_dates = sort_data(glob.glob(indir + "*aia*" + globexp + ".fits"))

    # find datetimes that are in *all* lists
    common_dates = list(set.intersection(*map(set, [con_dates, mag_dates, dop_dates, aia_dates])))

    # remove epochs that are missing in any data set from all data sets
    con_files = [con_files[idx] for idx, date in enumerate(con_dates) if date in common_dates]
    mag_files = [mag_files[idx] for idx, date in enumerate(mag_dates) if date in common_dates]
    dop_files = [dop_files[idx] for idx, date in enumerate(dop_dates) if date in common_dates]
    aia_files = [aia_files[idx] for idx, date in enumerate(aia_dates) if date in common_dates]
    return con_files, mag_files, dop_files, aia_files

def sort_data(f_list):
    # sort, and only take unique dates
    dates, inds = np.unique(get_dates(f_list), return_index=True)
    return [f_list[i] for i in inds], dates

def get_date(f):
    if "aia" in f:
        s = re.search(r'\d{4}_\d{2}_\d{2}t\d{2}_\d{2}_\d{2}', f)
    else:
        s = re.search(r'\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}', f)

    # standardize time formats
    s = s.group()
    s = s.replace("t", "_")
    return round_time(date=dt.datetime.strptime(s, "%Y_%m_%d_%H_%M_%S"))

def get_dates(files):
    return list(map(get_date, files))

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

def organize_IO(indir, datadir=None, clobber=False, globexp=""):
    # find the input data and check the lengths
    assert isdir(indir)
    con_files, mag_files, dop_files, aia_files = find_data(indir, globexp=globexp)
    assert (len(con_files) == len(mag_files) == len(dop_files) == len(aia_files))

    # figure out data directories
    if datadir == None:
        globdir = globexp.replace("*","")
        datadir = str(root / "data") + "/" + globdir + "/"

    if not isdir(datadir):
        os.mkdir(datadir)

    # name output files
    fname1 = datadir + "thresholds.csv"
    fname2 = datadir + "pixel_stats.csv"
    fname3 = datadir + "light_stats.csv"
    fname4 = datadir + "velocities.csv"
    fname5 = datadir + "mag_stats.csv"
    fname6 = datadir + "unweighted_velocities.csv"
    fname7 = datadir + "intensities.csv"
    fname8 = datadir + "average_velocities.csv"

    # headers for output files
    header1 = ["mjd", "aia_thresh", "a_aia", "b_aia", "c_aia", "hmi_thresh1", "hmi_thresh2", "a_hmi", "b_hmi", "c_hmi"]
    header2 = ["mjd", "lo_mu", "hi_mu", "ffactor", "umb_frac", "blu_pen_frac", "red_pen_frac", "quiet_frac", "network_frac", "plage_frac"]
    header3 = ["mjd", "lo_mu", "hi_mu", "umb_frac", "blu_pen_frac", "red_pen_frac", "quiet_frac", "network_frac", "plage_frac"]
    header4 = ["mjd", "region", "lo_mu", "hi_mu", "v_hat", "v_phot", "v_quiet", "v_conv"]
    header5 = ["mjd", "region", "lo_mu", "hi_mu", "mag_avg", "mag_std", "mag_unsigned"]
    header6 = ["mjd", "region", "lo_mu", "hi_mu", "v_hat", "v_phot", "v_quiet", "v_conv"]
    header7 = ["mjd", "region", "lo_mu", "hi_mu", "avg_int", "avg_int_flat"]
    header8 = ["mjd", "region", "lo_mu", "hi_mu", "v_hat", "v_phot", "v_quiet", "v_conv"]

    # replace/create/modify output files
    fileset = (fname1, fname2, fname3, fname4, fname5, fname6, fname7, fname8)
    if clobber and any(map(exists, fileset)):
        # delete the files
        clean_output_directory(*fileset)

        # create the files with headers
        create_file(fname1, header1)
        create_file(fname2, header2)
        create_file(fname3, header3)
        create_file(fname4, header4)
        create_file(fname5, header5)
        create_file(fname6, header6)
        create_file(fname7, header7)
        create_file(fname8, header8)
    elif all(map(exists, fileset)) and all(map(lambda x: getsize(x) > 0, fileset)):
        # get list of dates from file
        mjd_list = find_all_dates(fname1)

        # convert to Time objects and round to nearest hour
        mjd_list = list(map(lambda x: Time(x, format="mjd"), mjd_list))
        mjd_list = list(map(lambda x: round_time(date=x.datetime), mjd_list))

        # subset the input data to list to only include dates not seen here
        common_dates = list(set.intersection(*map(set, [get_dates(con_files), mjd_list])))

        # remove epochs that are missing in any data set from all data sets
        con_files = [con_files[idx] for idx, date in enumerate(get_dates(con_files)) if date not in common_dates]
        mag_files = [mag_files[idx] for idx, date in enumerate(get_dates(mag_files)) if date not in common_dates]
        dop_files = [dop_files[idx] for idx, date in enumerate(get_dates(dop_files)) if date not in common_dates]
        aia_files = [aia_files[idx] for idx, date in enumerate(get_dates(aia_files)) if date not in common_dates]
    else:
        create_file(fname1, header1)
        create_file(fname2, header2)
        create_file(fname3, header3)
        create_file(fname4, header4)
        create_file(fname5, header5)
        create_file(fname6, header6)
        create_file(fname7, header7)
        create_file(fname8, header8)

    return con_files, mag_files, dop_files, aia_files

def clean_output_directory(*fnames):
    for fname in fnames:
        truncate_output_file(fname)
        fname_mp = glob.glob(split(fname)[0] + "/tmp/" + splitext(split(fname)[1])[0] + "_*")
        if not not fname_mp:
            for f_mp in fname_mp:
                os.remove(f_mp)
    return None

def truncate_output_file(*fnames):
    # truncate the file if it does exist
    for fname in fnames:
        if exists(fname):
            with open(fname, "w") as f:
                f.truncate()
    return None

def find_all_dates(fname):
    mjd_list = []
    with open(fname, "r") as f:
        for line in f:
            if "mjd" in line:
                continue
            mjd_list.append(line.split(",")[0])
    return mjd_list

def create_file(fname, header=None):
    with open(fname, "w") as f:
        writer = csv.writer(f)
        if header is not None:
            writer.writerow(header)
    return None

def write_results_to_file(fname, *args):
    assert exists(fname)

    # parse out args
    lines = [a for a in args]

    if any(isinstance(el, list) for el in lines):
        for el in lines:
            write_results_to_file(fname, *el)
    else:
        # write to disk
        with open(fname, "a") as f:
            writer = csv.writer(f)
            writer.writerow(lines)
    return None

def stitch_output_files(fname, files, delete=False):
    with open(fname, "a") as f:
        for file in files:
            with open(file, "r") as g:
                for line in g:
                    if "mjd" in line:
                        continue
                    f.write(line)

    if delete:
        for f in files:
            os.remove(f)
    return None
