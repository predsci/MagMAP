#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import datetime
import argparse
import pytz
import signal
#
import magmap.data.download.drms_helpers as drms_helpers
import magmap.utilities.file_io.io_helpers as io_helpers
#
########################################
#  MAGMAP_GET_DATA.PY:  Version 1.0.0  #
########################################
########################################################################
#          Predictive Science Inc.
#          www.predsci.com
#          San Diego, California, USA 92121
########################################################################
# Copyright 2024 Predictive Science Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied.
# See the License for the specific language governing permissions and
# limitations under the License.
########################################################################

def signal_handler(signal, frame):
  print('You pressed Ctrl+C! Stopping!')
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def argParsing():
  parser = argparse.ArgumentParser(description='magmap_get_data.py:  This script uses magmap to query/download magnetograms and writes them to disk.  It also generates an index csv file for the files.  It can run in "update" mode to only download the newly available files into a pre-exisiting magmap download and regenerate the index csv file.')

  parser.add_argument(
            'period_start_input',
            help='Start date for data aquisition in the format:  YYYY-MM-DDTHH:MM:SS')

  parser.add_argument(
            'period_end_input',
            help='End date for data aquisition in the format:  YYYY-MM-DDTHH:MM:SS')

  parser.add_argument('-odir',
            help='Full path of the directory for output. Default is ${PWD}/magmap_data_disks',
            dest='download_dir',
            default='magmap_data_disks',
            required=False)

  parser.add_argument('-cadence',
            help='Magnetogram download cadence (hours)',
            dest='input_interval_cadence',
            default=1.0,
            type=float,
            required=False)

  parser.add_argument('-search_width',
            help='Magnetogram download query window half-width (>0) (hours).',
            dest='input_del_interval',
            default=0.33,
            type=float,
            required=False)

  parser.add_argument('-index_filename',
            help='Filename for the CSV index file.  Default is "magmap_all.csv',
            dest='index_file',
            default='magmap_disks_all.csv',
            required=False)

  return parser.parse_args()

def run(args):

# ---- Inputs -----------------------------
  download_dir           = args.download_dir
  period_start_input     = args.period_start_input
  period_end_input       = args.period_end_input
  input_interval_cadence = args.input_interval_cadence
  input_del_interval     = args.input_del_interval
  index_file             = args.index_file
#
  if input_del_interval <= 0.0:
    print('WARNING!  Search width was set to <= 0.  Setting to 0.0167 hours.')
    input_del_interval = 0.0167
#
  if not os.path.exists(download_dir):
    print('Download directory does not exist, creating '+download_dir)
    os.makedirs(download_dir)
#
# ----- End Inputs -------------------------
#
# define image search interval cadence and width
  interval_cadence = datetime.timedelta(hours=input_interval_cadence)
  del_interval     = datetime.timedelta(hours=input_del_interval)

# select the data series [only 1 available for now]
  series = 'hmi.m_720s'
  filters = ['QUALITY=0', ]

# initial time
  period_start = datetime.datetime.strptime(period_start_input, '%Y-%m-%dT%H:%M:%S')
  period_start = pytz.utc.localize(period_start)

# ending time
  period_end = datetime.datetime.strptime(period_end_input, '%Y-%m-%dT%H:%M:%S')
  period_end = pytz.utc.localize(period_end)

# define target times over download period using interval_cadence (image times in astropy Time() format)
  target_times = pd.date_range(start=period_start, end=period_end, freq=interval_cadence).to_pydatetime()
  query_range = [period_start-del_interval, target_times[-1]+del_interval]

# initialize the helper class for HMI
  hmi = drms_helpers.HMI_M720s(series=series, filters=filters, verbose=True)

# query available magnetograms
  available_hmi = hmi.query_time_interval(time_range=query_range)

# generate DataFrame that defines synchronic target times as well as min/max limits
  match_times = pd.DataFrame({'target_time': target_times, 'hmi_time': target_times})
# match closest available to each target time
  available_datetimes = np.array([
    datetime.datetime.strptime(x, "%Y.%m.%d_%H:%M:%S").replace(tzinfo=pytz.UTC)
    for x in available_hmi.time])

  total_possible_files = len(match_times)

  for index, row in match_times.iterrows():
    time_diff = available_datetimes - row.target_time.to_pydatetime()
    time_diff = np.abs(time_diff)
    best_match = time_diff.argmin()
    if time_diff[best_match] <= del_interval:
      # download resulting magnetogram
      print(f'\nStep {index+1}/{total_possible_files}:  Acquiring data for {available_hmi.loc[best_match].time}, quality code: {available_hmi.loc[best_match].quality}')
      sub_dir, fname, exit_flag = hmi.download_image_fixed_format(
          data_series=available_hmi.loc[best_match], base_dir=download_dir,
          update=True, overwrite=False, verbose=True
      )
    else:
      print(f'\nStep {index}/{total_possible_files}: --> NO SUITABLE DATA FOR SELECTED INTERVAL AT: {row.hmi_time}')

# read the updated filesystem
  available_raw = io_helpers.read_db_dir(download_dir)

# write to csv
  sub_dir = os.path.join(download_dir, "index_files")
  index_full_path = os.path.join(sub_dir, index_file)
  if not os.path.isdir(sub_dir):
    os.makedirs(sub_dir)
  print(f'\nGenerating index CSV file: {index_full_path}')
  available_raw.to_csv(index_full_path, index=False,
                       date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')

  print(f'\nMagMAP data acquisition complete!\n')
def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()
