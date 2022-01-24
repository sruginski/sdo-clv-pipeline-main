# necessary modules
import sunpy as sp
import sunpy.map as spm

# read headers and data
def read_header(file):
    return sp.io.fits.header_to_fits(sp.io.fits.get_header(file)[1])

def read_data(file):
    return spm.Map(file).data

