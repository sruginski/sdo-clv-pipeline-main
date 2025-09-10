from sdo_clv_pipeline.sdo_download import *

# data_dir = os.path.abspath(os.path.join(os.getcwd(), "..", "data"))
data_dir = os.path.abspath("/mnt/ceph/users/mpalumbo/new_sdo_data")
print(data_dir)
download_data(series="720", email="mlp95@psu.edu", outdir=data_dir, 
              start="2019/01/01", end="2019/01/31", 
              sample=4, overwrite=False, progress=True)
