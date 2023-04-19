import glob
from sdo_pypline.sdo_download import download_data
from sdo_pypline.paths import root

def main():
    # define params
    series = "720"
    outdir = "/data/mlp95/sdo_720/"
    sample = 4
    overwrite = False
    progress = True

    # define start and end
    start12="2012/01/01"
    end12="2012/12/31"

    download_data(series=series, outdir=outdir, start=start12, end=end12,
                  sample=sample, overwrite=overwrite, progress=progress)

    # define start and end
    start13="2013/01/01"
    end13="2013/12/31"

    download_data(series=series, outdir=outdir, start=start13, end=end13,
                  sample=sample, overwrite=overwrite, progress=progress)

    # define start and end
    start14="2014/01/01"
    end14="2014/12/31"

    download_data(series=series, outdir=outdir, start=start14, end=end14,
                  sample=sample, overwrite=overwrite, progress=progress)

    # define start and end
    start15="2015/01/01"
    end15="2015/12/31"

    download_data(series=series, outdir=outdir, start=start15, end=end15,
                  sample=sample, overwrite=overwrite, progress=progress)


if __name__ == "__main__":
    main()
