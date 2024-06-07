"""
Determine last available map and last available raw magnetogram.
Then convert newer images to maps.
Update index file.
"""

import os
import numpy as np
import pandas as pd
import datetime
import time
from multiprocessing import Pool

import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.maps.hipft_prep as hipft_prep
import magmap.maps.util.map_manip as map_manip
import magmap.utilities.file_io.io_helpers as io_helpers
import magmap.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# select the data series
series = 'hmi.m_720s'
# series = 'hmi.m_720s_nrt'

if series == 'hmi.m_720s':
    raw_dirname = 'hmi_m720s'
    map_dirname = 'hmi_hipft'

elif series == 'hmi.m_720s_nrt':
    raw_dirname = 'hmi_m720s_nrt'
    map_dirname = 'hmi_hipft_m720s_nrt'

# data-file dirs
base_dir = f'/Volumes/extdata3/oft'
raw_data_dir = f"{base_dir}/raw_data/{raw_dirname}"
map_data_dir = f"{base_dir}/processed_maps/{map_dirname}"

# index file name
index_file = "all-files.csv"

# number of parallel processors and threads-per-processor
nprocs = 1
tpp = 5
print_timings = True

# High-res map grid specifications
map_nxcoord = 10240
map_nycoord = 5120
R0 = 1.

# reduced map grid specifications
reduced_nxcoord = 1024
reduced_nycoord = 512

# ----- End Inputs -------------------------

# read available magnetograms filesystem
print("\nReading filesystem from dir: " + raw_data_dir + "\n")
available_raw = io_helpers.read_db_dir(raw_data_dir)
if len(available_raw) == 0:
    raise BaseException(f'Could not find any files in {raw_data_dir}')

print("\nReading filesystem from dir: " + map_data_dir + "\n")
available_map = io_helpers.read_db_dir(map_data_dir)

# select all maps between these dates
if len(available_map) > 0:
    min_datetime_thresh = available_map.date.iloc[-1].to_pydatetime()
else:
    # if no maps have been processed, take the earliest time minus 10 mins
    min_datetime_thresh = available_raw.date.iloc[0].to_pydatetime() - datetime.timedelta(minutes=10)

# select only files/maps between thresh datetimes
keep_ind = (available_raw.date > min_datetime_thresh)
available_raw = available_raw.loc[keep_ind, :]


# initiate a pool of processors
# p_pool = Pool(nprocs)
p_pool = None

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

# setup reduced-map grid
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
    # rel_path = '2016/03/23/hmi_m_720s_20160323T155958.fits'
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

    # approximate 'full' map resolution
    # if index == 0:
    #     # image latitude per pixel at equator
    #     eq_radian_per_pixel = np.max(np.abs(np.diff(hmi_im.lat[1000:3000, 2049])))
    #     full_map_nycoord = np.ceil(np.pi/eq_radian_per_pixel)
    #     full_map_nxcoord = 2*full_map_nycoord

    # for debugging, plot the image and transformed CR coords
    # psi_plt.PlotDisk_wInterpGrid(disk_image=hmi_im, nfig=0, plot_attr="Br", map_x=x_axis, map_y=sin_lat, R0=R0)

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
# write to csv
index_full_path = os.path.join(map_data_dir, "index_files", index_file)
hipft_index.to_csv(index_full_path, index=False,
                   date_format="%Y-%m-%dT%H:%M:%S", float_format='%.5f')

# p_pool.close()

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

