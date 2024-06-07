"""
Specify times for HMI magnetogram download.
Query available images and download best matches.
"""

import os
import numpy as np
import pandas as pd
import datetime

import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.data.download.drms_helpers as drms_helpers
import magmap.maps.hipft_prep as hipft_prep
import magmap.maps.util.map_manip as map_manip

# ---- Inputs -----------------------------

# Specify a vector of query times
period_start = datetime.datetime(year=2012, month=1, day=15, hour=0)
period_end = datetime.datetime(year=2012, month=1, day=17, hour=0)
# period_end = datetime.datetime(year=2012, month=1, day=2, hour=0)
period_range = [period_start, period_end]

# define image search interval cadence and width
# interval_cadence = datetime.timedelta(hours=1)
# del_interval = datetime.timedelta(minutes=20)
interval_cadence = datetime.timedelta(minutes=12)
del_interval = datetime.timedelta(minutes=6)
# define target times over download period using interval_cadence (image times in astropy Time() format)
target_times = pd.date_range(start=period_start, end=period_end, freq=interval_cadence).to_pydatetime()
query_range = [period_start-del_interval, period_end+del_interval]

# specify path and filename for download_results file
down_results_dir = "/Users/turtle/Dropbox/MyOFT/download_test/download_results"
download_results_filename = "download_results_" + period_start.strftime("%Y-%m-%dT%H_%M_%S") + ".csv"
results_path = os.path.join(down_results_dir, download_results_filename)
hipft_text_filename = "hipft_input_" + period_start.strftime("%Y-%m-%dT%H_%M_%S") + ".csv"
hipft_text_path = os.path.join(down_results_dir, hipft_text_filename)

# data-file dirs
raw_data_dir = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw"
# map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"

# ----- End Inputs -------------------------

# initialize the helper class for HMI
hmi = drms_helpers.HMI_Mrmap_latlon_720s(verbose=True)

# query available magnetograms
available_hmi = hmi.query_time_interval(time_range=query_range)

# generate DataFrame that defines synchronic target times as well as min/max limits
match_times = pd.DataFrame({'target_time': target_times, 'hmi_time': target_times,
                            'url': ["", ]*len(target_times), 'raw_path': None,
                            'map_path':None
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



