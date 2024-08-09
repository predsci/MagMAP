
############################################################
# Zip all OFTpy maps by year
# and include an index file
#
#
#
############################################################

import os
import pandas as pd
import shutil

# data base directory
root_dir = "/Volumes/extdata3/oft/processed_maps/hmi_hipft"
# zip dir
zip_dir = "/Volumes/extdata3/oft/processed_maps/hmi_hipft_zips"
# set year to be compressed
year = 2022


# load index file
all_maps = pd.read_csv("/Volumes/extdata3/oft/processed_maps/hmi_hipft/index_files/all-files.csv")
# convert to timestamp
# all_maps.obs_datetime_utc = pd.to_datetime(all_maps.obs_datetime_utc, format='%Y-%m-%dT%H:%M:%S')
all_maps.obs_datetime_utc = pd.to_datetime(all_maps.obs_datetime_utc, format='%Y-%m-%d %H:%M:%S')
# subset index file by year
all_maps_year = all_maps.obs_datetime_utc.dt.year
all_maps = all_maps.loc[all_maps_year == year, :]
# write the single-year index file as a csv
year_index_file = os.path.join(zip_dir, str(year) + "_file_index.csv")
all_maps.to_csv(year_index_file, index=False, date_format="%Y-%m-%d %H:%M:%S", float_format='%.5f')

# set archive directory
archive_dir = str(year)
# set zipname
zip_fname = str(year)
zip_path = os.path.join(zip_dir, zip_fname)
print("Zipping", archive_dir, "to", zip_path)
shutil.make_archive(zip_path, 'gztar', root_dir=root_dir, base_dir=archive_dir)







