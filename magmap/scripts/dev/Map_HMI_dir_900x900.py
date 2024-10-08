"""
Specify times for HMI magnetogram download.
Query available images and download best matches.
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
# data-file dirs
raw_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_raw"
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800"
map_cor_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_cor"
helio_proj_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj"
helio_proj_post_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_postBr"
helio_proj_noBr_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_noBr"
helio_proj_sqrtBr_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_sqrtBr"
# raw_data_dir = "/Users/turtle/data/oft/hmi_raw"
# map_data_dir = "/Users/turtle/data/oft/hmi_map"

# select all maps between these dates
min_datetime_thresh = datetime.datetime(year=2012, month=1, day=15, hour=0)
max_datetime_thresh = datetime.datetime(year=2012, month=1, day=16, hour=0)

# number of processors for interpolation parallelization
nprocs = 4
tpp = 5

# High-res map grid specifications
map_nxcoord = 10240
map_nycoord = 5120
R0 = 1.

# reduced map grid specifications
reduced_nxcoord = 1800
reduced_nycoord = 900

# ----- End Inputs -------------------------

# initiate a pool of processors
p_pool = Pool(nprocs)

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

# setup reduced-map grid (specific to hmi_Mrmap_latlon_720s grid),
# with added points at poles and periodic boundary to match OFTpy
# mapping standard
reduced_x = np.linspace(0.1, 359.9, reduced_nxcoord) * np.pi/180
reduced_x = np.append(0., reduced_x)
reduced_x = np.append(reduced_x, 2*np.pi)
reduced_y = np.linspace(-89.9, 89.9, reduced_nycoord) * np.pi/180
reduced_y = np.append(-np.pi/2, reduced_y)
reduced_y = np.append(reduced_y, np.pi/2)
# interp expects sin(lat)
reduced_sin_lat = np.sin(reduced_y)

# read available magnetograms filesystem
print("\nReading filesystem from dir: " + raw_data_dir + "\n")

available_raw = io_helpers.read_db_dir(raw_data_dir)

# select only files/maps between thresh datetimes
keep_ind = (available_raw.date >= min_datetime_thresh) & \
           (available_raw.date <= max_datetime_thresh)
available_raw = available_raw.loc[keep_ind, :]


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
    # # for the purpose of this script, skip if file already exists
    # if os.path.exists(os.path.join(map_data_dir, map_rel)):
    #     print("Map file already exists. SKIPPING!")
    #     continue
    # load to LosMagneto object
    full_path = os.path.join(raw_data_dir, rel_path)
    hmi_im = psi_dtypes.read_hmi720s(full_path, make_map=False, solar_north_up=False)
    IOtime += time.time() - start_time

    # convert to Br (with slightly larger radius to pad image for interpolation)
    start_time = time.time()
    hmi_im.add_Br(mu_thresh=0.01, R0=R0)
    image_proc_time += time.time() - start_time

    start_time = time.time()
    # interpolate to map
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="Br",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool)
    interp_time += time.time() - start_time

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

    ## --- repeat with radius correction (jitter test) -----------------------
    if not os.path.exists(os.path.join(map_cor_dir, sub_dir)):
        os.makedirs(os.path.join(map_cor_dir, sub_dir), mode=0o755)
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=True)

    # convert interpolated map values to Br
    data_index = hmi_map.data > hmi_map.no_data_val
    hmi_map.data[data_index] = hmi_map.data[data_index] / hmi_map.mu[data_index]

    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                              periodic_x=True, y_units='sinlat', uniform_poles=True,
                                              uniform_no_data=True)

    # assign map y-axis back to phi
    reduced_map.y = reduced_y
    # set assimilation weights
    reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

    # write to hipft file
    reduced_map.write_to_file(map_cor_dir, map_type='magneto', filename=map_rel)

    ## --- repeat with radius correction (helio-projective) -----------------------
    if not os.path.exists(os.path.join(helio_proj_dir, sub_dir)):
        os.makedirs(os.path.join(helio_proj_dir, sub_dir), mode=0o755)
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="Br",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=False, helio_proj=True)

    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                              periodic_x=True, y_units='sinlat', uniform_poles=True,
                                              uniform_no_data=True)

    # assign map y-axis back to phi
    reduced_map.y = reduced_y
    # set assimilation weights
    reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

    # write to hipft file
    reduced_map.write_to_file(helio_proj_dir, map_type='magneto', filename=map_rel)

    ## --- repeat with radius correction (helio-projective) and post Br calc -----------------------
    if not os.path.exists(os.path.join(helio_proj_post_dir, sub_dir)):
        os.makedirs(os.path.join(helio_proj_post_dir, sub_dir), mode=0o755)
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=False, helio_proj=True)

    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                              periodic_x=True, y_units='sinlat', uniform_poles=True,
                                              uniform_no_data=True)

    # assign map y-axis back to phi
    reduced_map.y = reduced_y
    # convert interpolated map values to Br
    data_index = reduced_map.data > reduced_map.no_data_val
    reduced_map.data[data_index] = reduced_map.data[data_index] / reduced_map.mu[data_index]
    # set assimilation weights
    reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

    # write to hipft file
    reduced_map.write_to_file(helio_proj_post_dir, map_type='magneto', filename=map_rel)

    print("")

    ## --- repeat with radius correction (helio-projective) and no Br calc -----------------------
    if not os.path.exists(os.path.join(helio_proj_noBr_dir, sub_dir)):
        os.makedirs(os.path.join(helio_proj_noBr_dir, sub_dir), mode=0o755)
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=False, helio_proj=True)

    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                              periodic_x=True, y_units='sinlat', uniform_poles=True,
                                              uniform_no_data=True)

    # assign map y-axis back to phi
    reduced_map.y = reduced_y

    # set assimilation weights
    reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

    # write to hipft file
    reduced_map.write_to_file(helio_proj_noBr_dir, map_type='magneto', filename=map_rel)

    print("")

    ## --- repeat with radius correction (helio-projective) and sqrt(mu) Br calc -----------------------
    if not os.path.exists(os.path.join(helio_proj_sqrtBr_dir, sub_dir)):
        os.makedirs(os.path.join(helio_proj_sqrtBr_dir, sub_dir), mode=0o755)
    hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                   nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=False, helio_proj=True)

    # down-sample by integration
    reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                              periodic_x=True, y_units='sinlat', uniform_poles=True,
                                              uniform_no_data=True)

    # assign map y-axis back to phi
    reduced_map.y = reduced_y
    # convert interpolated map values to Br using 1/sqrt(mu)
    data_index = (reduced_map.data > reduced_map.no_data_val) & \
                 (reduced_map.mu >= 0.)
    reduced_map.data[data_index] = reduced_map.data[data_index] / np.sqrt(reduced_map.mu[data_index])

    # set assimilation weights
    reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

    # write to hipft file
    reduced_map.write_to_file(helio_proj_sqrtBr_dir, map_type='magneto', filename=map_rel)

    print("")

p_pool.close()

n_its = available_raw.shape[0]
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

