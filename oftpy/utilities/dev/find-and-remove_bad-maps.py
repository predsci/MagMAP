
import datetime
import pandas as pd
import numpy as np
import os

import oftpy.utilities.file_io.io_helpers as io_helpers
import oftpy.utilities.datatypes.datatypes as psi_dt

# ---- Inputs -----------------------------
# data-file dirs
raw_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_raw"
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"
del_summary_path = "/Volumes/terminus_ext/HMI_M720s"


# ---- Main -------------------------------
print("\nReading filesystem from dir: " + map_data_dir + "\n")
# get a list of raw filenames
available_maps = io_helpers.read_db_dir(map_data_dir)

deleted_index = list()
for index, row in available_maps.iterrows():
    # open the map
    full_path = os.path.join(map_data_dir, row.rel_path)
    if (index % 100) == 0:
        print("Opening map: ", full_path, sep="")
    check_map = psi_dt.read_hipft_map(full_path)
    # check for nans in 'data' attribute
    if np.isnan(check_map.data).any():
        # if there are nans, delete the map
        os.remove(full_path)
        deleted_index.append(index)
        print("\nDeleting bad file: ", full_path, "\n", sep="")

# reduce map list to only those deleted
del_maps = available_maps.loc[deleted_index]
cur_time = datetime.datetime.today()
del_filename = "MapDelSummary_" + cur_time.strftime("%Y-%m-%d_%H_%M") + ".csv"

del_file_path = os.path.join(del_summary_path, del_filename)
del_maps.to_csv(del_file_path, index=False,
                date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')
