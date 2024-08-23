"""
Read existing filesystem.
Query/Download available images from latest local file to current timestamp.
Re-read filesystem and generate index csv

NOTES on hmi NRT magnetograms
  - NOTE for NRT it seems that 1024 is applied to ALL files --> probably the bit for NRT
    HOWEVER, comparing the quality over time to science data on 2024/01/21 and 2024/01/22
    it seems that 72704 in NRT will correspond to 0 in Science and they also look similar
    --> for now include both using the OR statement
- i was able to examine this with the following settings:
    interval_cadence = datetime.timedelta(minutes=12)
    del_interval = datetime.timedelta(minutes=6)
    period_start = datetime.datetime(2024, 1, 22, 1, 0, 0, 1, tzinfo=datetime.timezone.utc)
    period_end = datetime.datetime(2024, 1, 22, 2, 0, 0, 1, tzinfo=datetime.timezone.utc)
"""

import os
import numpy as np
import pandas as pd
import datetime

import pytz

import magmap.data.download.drms_helpers as drms_helpers
import magmap.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------
# select the data series (see notes at the top about the filter)
series = 'hmi.m_720s'
# series = 'hmi.m_720s_nrt'

if series == 'hmi.m_720s':
    raw_dirname = 'hmi_m720s'
    filters = ['QUALITY=0', ]

elif series == 'hmi.m_720s_nrt':
    raw_dirname = 'hmi_m720s_nrt'
    filters = ['QUALITY=1024 or QUALITY=72704', ]

# define image search interval cadence and width
interval_cadence = datetime.timedelta(hours=1)
del_interval = datetime.timedelta(minutes=20)

# index file name
index_file = "all-files.csv"

# data-file dirs
base_dir = f'/Volumes/extdata3/oft'
raw_data_dir = f"{base_dir}/raw_data/{raw_dirname}"

# ----- End Inputs -------------------------

# read data directory
available_raw = io_helpers.read_db_dir(raw_data_dir)

# define date range to search

# initial time
if len(available_raw) > 0:
    period_start = (available_raw.date.iloc[-1].to_pydatetime()).astimezone(datetime.timezone.utc)
else:
    period_start = datetime.datetime(2024, 1, 1, 0, 0, 0, 1, tzinfo=datetime.timezone.utc)

# optionally override the start time to force a re-examination
# period_start = datetime.datetime(2024, 1, 1, 0, 0, 0, 1, tzinfo=datetime.timezone.utc)

# ending time
period_end = datetime.datetime.now(datetime.timezone.utc)

period_range = [period_start, period_end]

# define target times over download period using interval_cadence (image times in astropy Time() format)
rounded_start = period_start.replace(second=0, microsecond=0, minute=0, hour=period_start.hour) + \
                datetime.timedelta(hours=period_start.minute//30)
target_times = pd.date_range(start=rounded_start, end=period_end, freq=interval_cadence).to_pydatetime()
query_range = [period_start-del_interval, period_end]

# initialize the helper class for HMI
hmi = drms_helpers.HMI_M720s(verbose=True, series=series, filters=filters)

# query available magnetograms
available_hmi = hmi.query_time_interval(time_range=query_range)

# generate DataFrame that defines synchronic target times as well as min/max limits
match_times = pd.DataFrame({'target_time': target_times, 'hmi_time': target_times})
# match closest available to each target time
available_datetimes = np.array([
    datetime.datetime.strptime(x, "%Y.%m.%d_%H:%M:%S").replace(tzinfo=pytz.UTC)
    for x in available_hmi.time])

for index, row in match_times.iterrows():
    time_diff = available_datetimes - row.target_time.to_pydatetime()
    time_diff = np.abs(time_diff)
    best_match = time_diff.argmin()
    if time_diff[best_match] <= del_interval:
        # download resulting magnetogram
        print(f'Acquiring data for time: {available_hmi.loc[best_match].time}, quality: {available_hmi.loc[best_match].quality}')
        sub_dir, fname, exit_flag = hmi.download_image_fixed_format(
            data_series=available_hmi.loc[best_match], base_dir=raw_data_dir,
            update=True, overwrite=False, verbose=True
        )
    else:
        print(f'  NO SUITABLE DATA FOR INTERVAL AT: {row.hmi_time}')

# read the updated filesystem
available_raw = io_helpers.read_db_dir(raw_data_dir)

# write to csv
index_full_path = os.path.join(raw_data_dir, "index_files", index_file)
available_raw.to_csv(index_full_path, index=False,
                     date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')
