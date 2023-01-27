"""
Specify times for HMI magnetogram download.
Query available images and download best matches.
"""

import os
import numpy as np
import pandas as pd
import datetime

import oftpy.utilities.datatypes.datatypes as psi_dtypes
import oftpy.data.download.drms_helpers as drms_helpers
import oftpy.maps.hipft_prep as hipft_prep
import oftpy.maps.util.map_manip as map_manip

# ---- Inputs -----------------------------

# Specify a vector of query times
# Specify a vector of query times
period_start = datetime.datetime(year=2010, month=1, day=1, hour=0)
period_end = datetime.datetime(year=2011, month=1, day=1, hour=0)
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
# raw_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_raw"
raw_data_dir = "/Volumes/extdata3/oft/raw_data/hmi_m720s"
# map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"

# High-res map grid specifications
map_nxcoord = 10240
map_nycoord = 5120
R0 = 1.

# reduced map grid specifications
reduced_nxcoord = 1024
reduced_nycoord = 512

# ----- End Inputs -------------------------

# setup map grid
y_range = [-np.pi/2, np.pi/2]
x_range = [0, 2 * np.pi]
x_axis = np.linspace(x_range[0], x_range[1], map_nxcoord)
y_axis = np.linspace(y_range[0], y_range[1], map_nycoord)
# interp expects sin(lat)
sin_lat = np.sin(y_axis)

# setup reduced-map grid
reduced_x = np.linspace(x_range[0], x_range[1], reduced_nxcoord)
reduced_y = np.linspace(y_range[0], y_range[1], reduced_nycoord)
# interp expects sin(lat)
reduced_sin_lat = np.sin(reduced_y)

# initialize the helper class for HMI
hmi = drms_helpers.HMI_M720s(verbose=True)

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
#
# # write download results to file
# match_times.to_csv(results_path, index=False, date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')
#
# ## FORMAT SUMMARY FILE FOR HIPFT
# # remove rows with hmi_time==None
# reduced_match_times = match_times.loc[~match_times.hmi_time.isna(), :]
# reduced_match_times = reduced_match_times.reset_index()
# # save summary dataframe to file
# write_df = reduced_match_times.loc[:, ['target_time', 'hmi_time', 'map_path']]
# write_df = write_df.rename(columns=dict(target_time='target_datetime', hmi_time='hmi_datetime'))
# # add fractional days since unix-epoch
# target_datetime = write_df.target_datetime.dt.to_pydatetime()
# target_unix_seconds = [float(target_datetime[ii].strftime("%s")) for ii in range(len(target_datetime))]
# target_unix_days = [x/(60*60*24) for x in target_unix_seconds]
# hmi_datetime = write_df.hmi_datetime.dt.to_pydatetime()
# hmi_unix_seconds = [float(hmi_datetime[ii].strftime("%s")) for ii in range(len(hmi_datetime))]
# hmi_unix_days = [x/(60*60*24) for x in hmi_unix_seconds]
# # add new columns to dataframe
# unix_time_df = pd.DataFrame(dict(target_unix_days=target_unix_days, hmi_unix_days=hmi_unix_days))
# write_df = pd.concat([unix_time_df, write_df], axis=1, ignore_index=True)
# # write to csv
# write_df.to_csv(hipft_text_path, index=False, date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')


