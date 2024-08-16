"""
Read existing filesystem.
Query/Download available images from latest local file to current timestamp.
Re-read filesystem and generate index csv


"""

import os
import numpy as np
import pandas as pd
import datetime

import pytz

import magmap.data.download.drms_helpers as drms_helpers
import magmap.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------
# --- required ---
download_dir = '/Users/turtle/data/oft/test/hmi_m720s'   # This script will create YYYY/MM/DD/ sub-directories within download_dir
period_start_input = "2024/01/01T00:00:00"
period_end_input = "2024/01/02T00:00:00"

# --- optional ---
input_interval_cadence = 1    # Image download cadence (hours)
input_del_interval = 20       # Query window half-width (minutes). Generally should be >0, and < 0.5*interval_cadence

# index file name
index_file = "all-files.csv"  # The name of the index file. Will be saved to download_dir + /index_files/ + index_file

# ----- End Inputs -------------------------

# define image search interval cadence and width
interval_cadence = datetime.timedelta(hours=input_interval_cadence)
del_interval = datetime.timedelta(minutes=input_del_interval)

# select the data series (see notes at the top about the filter)
series = 'hmi.m_720s'
filters = ['QUALITY=0', ]

# initial time
period_start = datetime.datetime.strptime(period_start_input, '%Y/%m/%dT%H:%M:%S')
period_start = pytz.utc.localize(period_start)

# ending time
period_end = datetime.datetime.strptime(period_end_input, '%Y/%m/%dT%H:%M:%S')
period_end = pytz.utc.localize(period_end)

period_range = [period_start, period_end]

# define target times over download period using interval_cadence (image times in astropy Time() format)
# rounded_start = period_start.replace(second=0, microsecond=0, minute=0, hour=period_start.hour) + \
#                 datetime.timedelta(hours=period_start.minute//30)
# target_times = pd.date_range(start=rounded_start, end=period_end, freq=interval_cadence).to_pydatetime()
target_times = pd.date_range(start=period_start, end=period_end, freq=interval_cadence).to_pydatetime()
query_range = [period_start-del_interval, target_times[-1]+del_interval]

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
            data_series=available_hmi.loc[best_match], base_dir=download_dir,
            update=True, overwrite=False, verbose=True
        )
    else:
        print(f'  NO SUITABLE DATA FOR INTERVAL AT: {row.hmi_time}')

# read the updated filesystem
available_raw = io_helpers.read_db_dir(download_dir)

# write to csv
sub_dir = os.path.join(download_dir, "index_files")
index_full_path = os.path.join(sub_dir, index_file)
if not os.path.isdir(sub_dir):
    os.makedirs(sub_dir)
available_raw.to_csv(index_full_path, index=False,
                     date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')


