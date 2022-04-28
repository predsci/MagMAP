
import pandas as pd
import datetime
import os

import oftpy.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------
# data-file dirs
raw_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_raw"
# map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"
# hipft_text_path = "/Volumes/terminus_ext/HMI_M720s"
map_data_dir = "/Volumes/terminus_ext/HMI_ron/hmi_map"
hipft_text_path = "/Volumes/terminus_ext/HMI_ron"


# select all maps between these dates
min_datetime_thresh = datetime.datetime(2010, 1, 1, 0, 0, 0)
max_datetime_thresh = datetime.datetime(2016, 2, 1, 0, 0, 0)

# ---- Main -------------------------------
# get a list of raw filenames
available_maps = io_helpers.read_db_dir(map_data_dir)


# save summary dataframe to file
write_df = available_maps.rename(columns=dict(date='hmi_datetime', rel_path='map_path'))
# write_df = write_df.loc[write_df.hmi_datetime < min_datetime_thresh, :]
keep_ind = (write_df.hmi_datetime >= min_datetime_thresh) & \
           (write_df.hmi_datetime <= max_datetime_thresh)
write_df = write_df.loc[keep_ind, :]
write_df.loc[:, 'target_datetime'] = write_df.hmi_datetime.dt.round('H')
# re-order columns and reset index
write_df = write_df.loc[:, ['target_datetime', 'hmi_datetime', 'map_path']]
write_df.reset_index(drop=True, inplace=True)

# add fractional days since unix-epoch
target_datetime = write_df.target_datetime.dt.to_pydatetime()
target_unix_seconds = [float(target_datetime[ii].strftime("%s")) for ii in range(len(target_datetime))]
target_unix_days = [x/(60*60*24) for x in target_unix_seconds]
hmi_datetime = write_df.hmi_datetime.dt.to_pydatetime()
hmi_unix_seconds = [float(hmi_datetime[ii].strftime("%s")) for ii in range(len(hmi_datetime))]
hmi_unix_days = [x/(60*60*24) for x in hmi_unix_seconds]
# add new columns to dataframe
unix_time_df = pd.DataFrame(dict(target_unix_days=target_unix_days, hmi_unix_days=hmi_unix_days))
write_df = pd.concat([unix_time_df, write_df], axis=1)

# generate a filename
map_index_filename = "maps-up-to_" + write_df.hmi_datetime.max().strftime("%Y-%m-%d") + ".csv"

# write to csv
write_df.to_csv(os.path.join(hipft_text_path, map_index_filename), index=False,
                date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')



