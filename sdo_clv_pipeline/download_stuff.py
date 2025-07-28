from sdo_clv_pipeline.sdo_download import *

data_dir = os.path.abspath(os.path.join(os.getcwd(), "..", "data"))
print(data_dir)
download_data(series="720", email="mlp95@psu.edu", outdir=data_dir, start="2014/02/01", end="2014/03/31", 
              sample=24, overwrite=False, progress=False)
