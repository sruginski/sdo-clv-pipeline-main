# necessary modules
import sunpy as sp
from astropy.io import fits

# read headers and data
def read_header(file):
    return fits.getheader(file, 1)

def read_data(file):
    return fits.getdata(file, 1)

