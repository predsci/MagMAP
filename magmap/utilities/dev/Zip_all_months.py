
############################################################
# Zip all OFTpy maps by month
#
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

# load index file
all_maps = pd.read_csv("/Volumes/extdata3/oft/processed_maps/hmi_hipft/index_files/all-files.csv")
# convert to timestamp
all_maps.obs_datetime_utc = pd.to_datetime(all_maps.obs_datetime_utc, format='%Y-%m-%dT%H:%M:%S')
all_maps['year'] = all_maps.obs_datetime_utc.dt.year
all_maps['month'] = all_maps.obs_datetime_utc.dt.month
# get all year/month combinations
all_months = all_maps.groupby(['year', 'month']).size().reset_index(name='count')

for index, row in all_months.iterrows():
    # archive_dir = os.path.join(root_dir, str(row.year), f"{row.month:02d}")
    archive_dir = os.path.join(str(row.year), f"{row.month:02d}")
    zip_fname = str(row.year) + "_" + f"{row.month:02d}"
    zip_path = os.path.join(zip_dir, zip_fname)
    print("Zipping", archive_dir, "to", zip_path)
    shutil.make_archive(zip_path, 'zip', root_dir=root_dir, base_dir=archive_dir)







