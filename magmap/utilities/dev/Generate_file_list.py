
# Read the file tree from designated directory and write a csv summary for HIPFT.

import pandas as pd
import datetime
import os
import astropy.time as astro_time

import magmap.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------
# data-file dirs
# raw_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_raw"
map_data_dir = "/Volumes/extdata3/oft/processed_maps/hmi_hipft"
hipft_text_path = "/Volumes/extdata3/oft/processed_maps/hmi_hipft/index_files"
# map_data_dir = "/Volumes/terminus_ext/HMI_ron/hmi_map"
# hipft_text_path = "/Volumes/terminus_ext/HMI_ron"


# select all maps between these dates
min_datetime_thresh = datetime.datetime(2010, 12, 31, 23, 31, 0, tzinfo=datetime.timezone.utc)
max_datetime_thresh = datetime.datetime(2022, 1, 1, 0, 0, 0, tzinfo=datetime.timezone.utc)
# or ignore date range and include all available files
all_flag = True

# ---- Main -------------------------------
# get a list of raw filenames
available_maps = io_helpers.read_db_dir(map_data_dir)


# save summary dataframe to file
write_df = available_maps.rename(columns=dict(date='obs_datetime_utc', rel_path='map_path'))
if not all_flag:
    keep_ind = (write_df.obs_datetime_utc >= min_datetime_thresh) & \
               (write_df.obs_datetime_utc <= max_datetime_thresh)
    write_df = write_df.loc[keep_ind, :]
write_df.loc[:, 'target_datetime_utc'] = write_df.obs_datetime_utc.dt.round('H')
# re-order columns and reset index
write_df = write_df.loc[:, ['target_datetime_utc', 'obs_datetime_utc', 'map_path']]
write_df.reset_index(drop=True, inplace=True)

# add julian-days
obs_astro_time = astro_time.Time(write_df.obs_datetime_utc)
obs_jdays = obs_astro_time.jd
# add new columns to dataframe
jd_time_df = pd.DataFrame(dict(obs_jd=obs_jdays))
write_df = pd.concat([jd_time_df, write_df], axis=1)
# re-order columns
write_df = write_df.loc[:, ['target_datetime_utc', 'obs_datetime_utc', 'obs_jd', 'map_path']]

# generate a filename
if all_flag:
    map_index_filename = "all-files.csv"
else:
    map_index_filename = "maps-up-to_" + write_df.obs_datetime_utc.max().strftime("%Y-%m-%d") + ".csv"

# write to csv
write_df.to_csv(os.path.join(hipft_text_path, map_index_filename), index=False,
                date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')



