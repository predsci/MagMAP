"""
Read existing filesystem.
Query/Download available images from latest local file to current timestamp.
Re-read filesystem and generate index csv
"""

import os
import numpy as np
import pandas as pd
import datetime

import oftpy.data.download.drms_helpers as drms_helpers
import oftpy.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------

# define image search interval cadence and width
interval_cadence = datetime.timedelta(hours=1)
del_interval = datetime.timedelta(minutes=20)

# index file name
index_file = "all-files.csv"

# data-file dirs
raw_data_dir = "/Volumes/extdata3/oft/raw_data/hmi_m720s"

# ----- End Inputs -------------------------

# read data directory
available_raw = io_helpers.read_db_dir(raw_data_dir)

# define date range to search
period_start = available_raw.date.iloc[-1].to_pydatetime()
period_end = datetime.datetime.now()
period_range = [period_start, period_end]

# define target times over download period using interval_cadence (image times in astropy Time() format)
target_times = pd.date_range(start=period_start, end=period_end, freq=interval_cadence).to_pydatetime()
query_range = [period_start-del_interval, period_end]

# initialize the helper class for HMI
hmi = drms_helpers.HMI_M720s(verbose=True)

# query available magnetograms
available_hmi = hmi.query_time_interval(time_range=query_range)

# generate DataFrame that defines synchronic target times as well as min/max limits
match_times = pd.DataFrame({'target_time': target_times, 'hmi_time': target_times,
                            'url': ["", ]*len(target_times), 'raw_path': None,
                            'map_path': None
                            })
# match closest available to each target time
available_datetimes = np.array([datetime.datetime.strptime(x, "%Y.%m.%d_%H:%M:%S") for x in available_hmi.time])

hmi_download = pd.DataFrame(columns=available_hmi.columns)

for index, row in match_times.iterrows():
    time_diff = available_datetimes - row.target_time.to_pydatetime()
    time_diff = np.abs(time_diff)
    best_match = time_diff.argmin()
    if time_diff[best_match] <= del_interval:
        # add to output dataframe
        match_times.loc[index, 'hmi_time'] = pd.Timestamp(available_datetimes[best_match])
        match_times.loc[index, 'url'] = available_hmi.url[best_match]
        hmi_download = hmi_download.append(available_hmi.loc[best_match])

        # download resulting magnetograms
        sub_dir, fname, exit_flag = hmi.download_image_fixed_format(
            data_series=available_hmi.loc[best_match], base_dir=raw_data_dir,
            update=True, overwrite=False, verbose=True
        )

    print("")

# read the updated filesystem
available_raw = io_helpers.read_db_dir(raw_data_dir)

# write to csv
index_full_path = os.path.join(raw_data_dir, "index_files", index_file)
available_raw.to_csv(index_full_path, index=False,
                     date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')
