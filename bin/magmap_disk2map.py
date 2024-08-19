"""
Convert LOS magnetogram disk images to B_r Carrington maps
Two modes:
  Automatic update - The script determines the last available map and
  last available raw magnetogram. All images that are newer than the
  most recent map are converted.
  Manual dates - (default) User inputs a start and end datetime.  All
  magnetograms in that range are converted.

"""

import os
import numpy as np
import pandas as pd
import datetime
import time
import pytz

import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.maps.hipft_prep as hipft_prep
import magmap.maps.util.map_manip as map_manip
import magmap.utilities.file_io.io_helpers as io_helpers

# ---- Inputs -----------------------------
# --- required ---
raw_data_dir = "/Users/turtle/data/oft/test/hmi_m720s"
map_data_dir = "/Users/turtle/data/oft/test/hmi_hipft"
period_start_input = "2024/01/01T00:00:00"  # To run in automatic mode, enter 'auto' for this variable
# period_start_input = "auto"
period_end_input = "2024/01/02T00:00:00"

# --- optional ---
# select the data series
series = 'hmi.m_720s'

# index file name
index_file = "all-files.csv"
# flag to turn timings on/off
print_timings = False

# High-res map grid specifications (interpolation grid)
map_nxcoord = 10240
map_nycoord = 5120
R0 = 1.

# reduced map grid specifications
reduced_nxcoord = 1024
reduced_nycoord = 512

# ----- End Inputs -------------------------

# read available magnetograms filesystem
if os.path.exists(os.path.join(raw_data_dir, "index_files", index_file)):
    available_raw = pd.read_csv(os.path.join(raw_data_dir, "index_files", index_file))
    available_raw['date'] = pd.to_datetime(available_raw['date'], format="%Y-%m-%dT%H:%M:%S").dt.tz_localize('UTC')
else:
    print("Did not find an index file at " + raw_data_dir + "/index_files/" + index_file)
    print("\nReading filesystem directly from dir: " + raw_data_dir + "\n")
    available_raw = io_helpers.read_db_dir(raw_data_dir)
if len(available_raw) == 0:
    raise BaseException(f'Could not find any files in {raw_data_dir}')

# initial time
if period_start_input.lower() == "auto":
    if os.path.exists(os.path.join(map_data_dir, "index_files", index_file)):
        available_map = pd.read_csv(os.path.join(map_data_dir, "index_files", index_file))
        available_map['date'] = pd.to_datetime(available_map['date'], format="%Y-%m-%dT%H:%M:%S").dt.tz_localize('UTC')
    else:
        print("\nReading filesystem directly from dir: " + map_data_dir + "\n")
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
    period_start = datetime.datetime.strptime(period_start_input, '%Y/%m/%dT%H:%M:%S')
    min_datetime_thresh = pytz.utc.localize(period_start)

    # ending time
    period_end = datetime.datetime.strptime(period_end_input, '%Y/%m/%dT%H:%M:%S')
    max_datetime_thresh = pytz.utc.localize(period_end)

    # select only raw files after the most recent map
    keep_ind = (available_raw.date >= min_datetime_thresh) & \
               (available_raw.date <= max_datetime_thresh)

# filter the available raw files to the selected set
available_raw = available_raw.loc[keep_ind, :]


# disable pool of processors option
p_pool = None
# number of parallel processors and threads-per-processor
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
sin_lat = np.sin(y_interp)

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


for index, row in available_raw.iterrows():
    start_time = time.time()
    rel_path = row.rel_path
    fname = os.path.basename(rel_path)
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
        print("Map file already exists. SKIPPING!")
        continue
    # load to LosMagneto object
    full_path = os.path.join(raw_data_dir, rel_path)
    hmi_im = psi_dtypes.read_hmi720s(full_path, make_map=False, solar_north_up=False)
    IOtime += time.time() - start_time

    start_time = time.time()
    # interpolate to map
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, no_data_val=-65500.,
                                   y_cor=False, helio_proj=True)
    interp_time += time.time() - start_time

    # convert interpolated map values to Br
    data_index = hmi_map.data > hmi_map.no_data_val
    hmi_map.data[data_index] = hmi_map.data[data_index] / hmi_map.mu[data_index]

    start_time = time.time()
    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
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
index_full_path = os.path.join(map_data_dir, "index_files", index_file)
hipft_index.to_csv(index_full_path, index=False,
                   date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')


n_its = available_raw.shape[0]
if n_its == 0:
    print("No maps to update.")
elif print_timings:
    IOtime = IOtime/n_its
    image_proc_time = image_proc_time/n_its
    interp_time = interp_time/n_its
    down_samp_time = down_samp_time/n_its
    map_proc_time = map_proc_time/n_its

    total_time = IOtime + image_proc_time + interp_time + down_samp_time + map_proc_time

    print("Total time: ", total_time)
    print("IO time: ", IOtime, "s or ", IOtime/total_time*100, "%", sep="")
    print("Image proc time: ", image_proc_time, "s or ", image_proc_time/total_time*100, "%", sep="")
    print("Interpolation time: ", interp_time, "s or ", interp_time/total_time*100, "%", sep="")
    print("Down sampling time: ", down_samp_time, "s or ", down_samp_time/total_time*100, "%", sep="")
    print("Map proc time: ", map_proc_time, "s or ", map_proc_time/total_time*100, "%", sep="")
