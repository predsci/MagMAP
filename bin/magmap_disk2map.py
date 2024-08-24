#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import datetime
import time
import argparse
import pytz
import signal
import warnings
#
import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.maps.hipft_prep as hipft_prep
import magmap.maps.util.map_manip as map_manip
import magmap.utilities.file_io.io_helpers as io_helpers
#
########################################
#  MAGMAP_DISK2MAP.PY:  Version 1.0.0  #
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
  parser = argparse.ArgumentParser(description='magmap_disk2map.py:  This script uses magmap to map magnetogram disk images (previously downloaded using magmap_get_data.py) onto a full-Sun map in Carrington coordinates.  It can run in "update" mode (the default mode) which will map all disk images listed in the chosen index file.  In this mode, if some disks are already mapped, it will only map the newly available disk images and regenerate the map index csv file.  In manual mode, a date range is entered and all disk images in that range are mapped.')

  parser.add_argument(
            'raw_data_dir',
            help='Directory of the previously downloaded disk magnetogram data set.')
             
  parser.add_argument('-odir',
            help='Full path of the directory for output. Default is ${PWD}/magmap_data_maps',
            dest='map_data_dir',
            default='magmap_data_maps',
            required=False)

  parser.add_argument('-startdate',
            help='Start date for disk mapping in the format:  YYYY-MM-DDTHH:MM:SS  If not specified, all disk images in the dataset+index will be mapped.',
            dest='period_start_input',
            required=False)
         
  parser.add_argument('-enddate',
            help='End date for disk mapping in the format:  YYYY-MM-DDTHH:MM:SS  If not specified, all disk images in the dataset+index will be mapped.',
            dest='period_end_input',
            required=False)

  parser.add_argument('-index_filename_disks',
            help='Filename for the CSV index file in the disk dataset.  Default is "magmap_disks_all.csv',
            dest='index_file_disks',
            default='magmap_disks_all.csv',
            required=False)

  parser.add_argument('-index_filename_maps',
            help='Filename for the CSV index file in the ouput map dataset.  Default is "magmap_maps_all.csv',
            dest='index_file_maps',
            default='magmap_maps_all.csv',
            required=False)

  parser.add_argument('-ntout',
            help='Colatitude grid size for output maps.',
            dest='reduced_nycoord',
            default=512,
            type=int,
            required=False)

  parser.add_argument('-npout',
            help='Longitude grid size for output maps.',
            dest='reduced_nxcoord',
            default=1024,
            type=int,
            required=False)

  parser.add_argument('-ntinterp',
            help='Colatitude grid size for high-res interpolation grid used in mapping.',
            dest='map_nycoord',
            default=5120,
            type=int,
            required=False)

  parser.add_argument('-npinterp',
            help='Longitude grid size for high-res interpolation grid used in mapping.',
            dest='map_nxcoord',
            default=10240,
            type=int,
            required=False)

  parser.add_argument('-r0',
            help='Assumed radius in units of solar radii (6.96e10cm).  Default is 1.0.',
            dest='R0',
            default=1.0,
            type=float,
            required=False)

#  parser.add_argument('-dataseries',
#            help='Data series of the downlaoded disk image dataset.',
#            dest='series',
#            default='hmi.m_720s',
#            type=float,
#            required=False)

  return parser.parse_args()

def run(args):

# ---- Inputs -----------------------------

  raw_data_dir       = args.raw_data_dir
  map_data_dir       = args.map_data_dir
  period_start_input = args.period_start_input
  period_end_input   = args.period_end_input
  index_file_disks   = args.index_file_disks
  index_file_maps    = args.index_file_maps
  reduced_nycoord    = args.reduced_nycoord
  reduced_nxcoord    = args.reduced_nxcoord
  map_nycoord        = args.map_nycoord
  map_nxcoord        = args.map_nxcoord
  R0                 = args.R0

  # Suppress Pandas FutureWarnings
  warnings.simplefilter(action='ignore', category=FutureWarning)

  if period_start_input is None and period_end_input is not None:
    print("Error!  You must either provide both a start and end data, OR do not specify either to use the update/auto mode.")
    return

  if period_start_input is not None and period_end_input is None:
    print("Error!  You must either provide both a start and end data, OR do not specify either to use the update/auto mode.")
    return

  if period_start_input is None and period_end_input is None:
     period_start_input = 'auto'

  if map_nycoord <= reduced_nycoord:
    print("Error!  The interpolation grid must be larger than the final grid.")
    return

  if map_nxcoord <= reduced_nxcoord:
    print("Error!  The interpolation grid must be larger than the final grid.")
    return
#
# ----- End Inputs -------------------------
#


# read available magnetograms filesystem
  if os.path.exists(os.path.join(raw_data_dir, "index_files", index_file_disks)):
      available_raw = pd.read_csv(os.path.join(raw_data_dir, "index_files", index_file_disks))
      available_raw['date'] = pd.to_datetime(available_raw['date'], format="%Y-%m-%dT%H:%M:%S").dt.tz_localize('UTC')
  else:
      print("Did not find an index file at " + raw_data_dir + "/index_files/" + index_file_disks)
      print("\nReading filesystem directly from dir: " + raw_data_dir)
      available_raw = io_helpers.read_db_dir(raw_data_dir)
  if len(available_raw) == 0:
      raise BaseException(f'Could not find any files in {raw_data_dir}')

# initial time
  if period_start_input.lower() == "auto":
    print("\nChecking if output map directory already exists...")
    if os.path.exists(os.path.join(map_data_dir, "index_files", index_file_maps)):
      available_map = pd.read_csv(os.path.join(map_data_dir, "index_files", index_file_maps))
      available_map['date'] = pd.to_datetime(available_map['target_datetime_utc'], format="%Y-%m-%dT%H:%M:%S").dt.tz_localize('UTC')
    else:
      print("\nNo index file found in output directory, constructing index info from filenames...")
      available_map = io_helpers.read_db_dir(map_data_dir)

    # select all maps between these dates
    if len(available_map) > 0:
      min_datetime_thresh = available_map.date.iloc[-1].to_pydatetime()
    else:
    # if no maps have been processed, take the earliest time minus 10 mins
      min_datetime_thresh = available_raw.date.iloc[0].to_pydatetime() - datetime.timedelta(minutes=10)

    # select only raw files after the most recent map
    keep_ind = (available_raw.date > min_datetime_thresh)

  else:
    # start time
    period_start = datetime.datetime.strptime(period_start_input, '%Y-%m-%dT%H:%M:%S')
    min_datetime_thresh = pytz.utc.localize(period_start)

    # ending time
    period_end = datetime.datetime.strptime(period_end_input, '%Y-%m-%dT%H:%M:%S')
    max_datetime_thresh = pytz.utc.localize(period_end)

    # select only raw files after the most recent map
    keep_ind = (available_raw.date >= min_datetime_thresh) & \
               (available_raw.date <= max_datetime_thresh)

# filter the available raw files to the selected set
  available_raw = available_raw.loc[keep_ind, :]

  print("\nDiscovered "+str(len(available_raw))+" disk images that need processing/mapping.")

# disable pool of processors option
  p_pool = None
# number of parallel processors and threads-per-processor
# this is not yet implemented.
  nprocs = 1
  tpp = 5

# setup map grid
  y_range = [-np.pi/2, np.pi/2]
  x_range = [0, 2 * np.pi]
  x_axis = np.linspace(x_range[0], x_range[1], map_nxcoord)
  y_axis = np.linspace(y_range[0], y_range[1], map_nycoord)
# adjust polar-coordinate interpolation-centers for half-pixel mesh
  dy = y_axis[1] - y_axis[0]
  y_interp = y_axis.copy()
  y_interp[0] = y_interp[0] + dy/4
  y_interp[-1] = y_interp[-1] - dy/4
# interp expects sin(lat)
  sin_lat_interp = np.sin(y_interp)
  sin_lat = np.sin(y_axis)

# setup reduced-map grid (for saving to file)
  reduced_x = np.linspace(x_range[0], x_range[1], reduced_nxcoord)
  reduced_y = np.linspace(y_range[0], y_range[1], reduced_nycoord)
# interp expects sin(lat)
  reduced_sin_lat = np.sin(reduced_y)

# initialize timing variables
  IOtime = 0
  image_proc_time = 0
  interp_time = 0
  down_samp_time = 0
  map_proc_time = 0
  loop_idx = 0

  for index, row in available_raw.iterrows():
 
    loop_idx = loop_idx + 1

    start_time = time.time()

    rel_path = row.rel_path

    fname   = os.path.basename(rel_path)
    sub_dir = os.path.dirname(rel_path)

    # determine path and filename
    map_filename = fname.replace("_m_", "_map_")
    map_filename = map_filename.replace(".fits", ".h5")
    map_rel = os.path.join(sub_dir, map_filename)
    # check that directory exists
    if not os.path.exists(os.path.join(map_data_dir, sub_dir)):
      os.makedirs(os.path.join(map_data_dir, sub_dir), mode=0o755)
    # for the purpose of this script, skip if file already exists
    if os.path.exists(os.path.join(map_data_dir, map_rel)):
      print("\nStep "+str(loop_idx)+"/"+str(len(available_raw))+":  Map file already exists. SKIPPING!")
      continue

    # load to LosMagneto object
    full_path = os.path.join(raw_data_dir, rel_path)

    hmi_im = psi_dtypes.read_hmi720s(full_path, make_map=False, solar_north_up=False)

    IOtime += time.time() - start_time

    start_time = time.time()

    print("\nStep "+str(loop_idx)+"/"+str(len(available_raw)))

    # interpolate to map
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat_interp, interp_field="data",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, no_data_val=-65500.,
                                   y_cor=False, helio_proj=True)
    # replace y-axis with original sin_lat grid
    hmi_map.y = sin_lat

    interp_time += time.time() - start_time

    # convert interpolated map values to Br
    data_index = hmi_map.data > hmi_map.no_data_val

    hmi_map.data[data_index] = hmi_map.data[data_index] / hmi_map.mu[data_index]

    start_time = time.time()
    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, 
                                              new_x=reduced_x, image_method=0,
                                              periodic_x=True, y_units='sinlat', uniform_poles=True,
                                              uniform_no_data=True)

    down_samp_time += time.time() - start_time

    start_time = time.time()
    # assign map y-axis back to phi
    reduced_map.y = reduced_y
    # set assimilation weights
    reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

    map_proc_time += time.time() - start_time

    start_time = time.time()
    # write to hipft file
    reduced_map.write_to_file(map_data_dir, map_type='magneto', filename=map_rel)

    IOtime += time.time() - start_time

    print("")

# read updated filesystem
  hipft_index = io_helpers.gen_hipft_index(map_data_dir)
# check that directory exists
  if not os.path.exists(os.path.join(map_data_dir, "index_files")):
    os.makedirs(os.path.join(map_data_dir, "index_files"), mode=0o755)
# write to csv
  index_full_path = os.path.join(map_data_dir, "index_files", index_file_maps)
  hipft_index.to_csv(index_full_path, index=False,
                     date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')


  n_its = available_raw.shape[0]
  if n_its > 0:
    IOtime = IOtime/n_its
    image_proc_time = image_proc_time/n_its
    interp_time = interp_time/n_its
    down_samp_time = down_samp_time/n_its
    map_proc_time = map_proc_time/n_its
    total_time = IOtime + image_proc_time + interp_time + down_samp_time + map_proc_time

    format_str = '%19s %10.2f sec  [%5.1f%%]'
    print(' ')
    print('----------------  Timing -------------------')
    print(' ')
    print(format_str % ("Input/Output:     ", IOtime,          IOtime/total_time*100))
    print(format_str % ("Image processing: ", image_proc_time, image_proc_time/total_time*100))
    print(format_str % ("Interpolation:    ", interp_time,     interp_time/total_time*100))
    print(format_str % ("Down sampling:    ", down_samp_time,  down_samp_time/total_time*100))
    print(format_str % ("Map processing:   ", map_proc_time,   map_proc_time/total_time*100))
    print('--------------------------------------------')
    print(format_str % ("Total:            ", total_time,      100.0))
    print('---------------------------------------------')
    
  print("")
  print("MagMAP mapping complete!")

def main():
  ## Get input agruments:
  args = argParsing()
  run(args)

if __name__ == '__main__':
  main()

