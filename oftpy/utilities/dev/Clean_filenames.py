
import os

import oftpy.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------
# data-file dirs
raw_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_raw"
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"

# ---- Main -------------------------------
# get a list of raw filenames
available_raw = io_helpers.read_db_dir(raw_data_dir)
for index, row in available_raw.iterrows():
    # apply filename correction
    old_rel_path = row.rel_path
    new_rel_path = old_rel_path.replace("_.", ".")
    # rename file
    os.rename(os.path.join(raw_data_dir, old_rel_path), os.path.join(raw_data_dir, new_rel_path))

available_maps = io_helpers.read_db_dir(map_data_dir)
for index, row in available_maps.iterrows():
    # apply filename correction
    old_rel_path = row.rel_path
    new_rel_path = old_rel_path.replace("_.", ".")
    # rename file
    os.rename(os.path.join(map_data_dir, old_rel_path), os.path.join(map_data_dir, new_rel_path))

