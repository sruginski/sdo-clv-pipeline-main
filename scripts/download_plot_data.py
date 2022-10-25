import glob
from sdo_pypline.sdo_download import download_data
from sdo_pypline.paths import root
datadir = str(root / "data") + "/"
plotdir = str(root / "figures") + "/"

def download_plot_data():
    # see if data are already downloaded
    files = glob.glob(datadir + "fits/*.fits")
    con = datadir + "fits/hmi_ic_45s_2014_01_07_00_01_30_tai_continuum.fits"
    mag = datadir + "fits/hmi_m_45s_2014_01_07_00_01_30_tai_magnetogram.fits"
    dop = datadir + "fits/hmi_v_45s_2014_01_07_00_01_30_tai_dopplergram.fits"
    aia = datadir + "fits/aia_lev1_1700a_2014_01_07t00_00_30_71z_image_lev1.fits"

    if any(map(lambda x: x not in files, (con, mag, dop, aia))):
        # download the data to plot
        print(">>> Downloading data")
        start = "2014/01/07"
        end = "2014/01/07"
        sample = 24
        con, mag, dop, aia = download_data(outdir=datadir + "fits/", start=start, end=end, sample=sample)

    return con, mag, dop, aia

if __name__ == "__main__":
    download_plot_data()
