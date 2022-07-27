
"""
Specify a single HMI magnetogram.
Convert to maps.  Save Hi-Res and HIPFT-Res maps.
"""

import os
import numpy as np
from multiprocessing import Pool

import oftpy.utilities.datatypes.datatypes as psi_dtypes
import oftpy.maps.hipft_prep as hipft_prep
import oftpy.maps.util.map_manip as map_manip
from oftpy.settings.info import DTypes


# ---- Inputs -----------------------------
fits_file_path = "/Volumes/extdata3/oft/raw_data/hmi_m720s/2014/06/15/hmi_m_720s_20140615T110008.fits"

# data-file dirs
raw_data_dir = "/Volumes/extdata3/oft/raw_data/hmi_m720s"
map_data_dir = "/Users/turtle/Dropbox/MyOFT/ron_tess"

# number of processors for interpolation parallelization
nprocs = 1
tpp = 5

# High-res map grid specifications
map_nxcoord = 10240
map_nycoord = 5120
R0 = 1.

# reduced map grid specifications
reduced_nxcoord = 1024
reduced_nycoord = 512

# ----- End Inputs -------------------------

# initiate a pool of processors
if nprocs > 1:
    p_pool = Pool(nprocs)
else:
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

# --- Main Script -------------------------------------------
rel_path = os.path.relpath(fits_file_path, raw_data_dir)
fname = os.path.basename(rel_path)
sub_dir = os.path.dirname(rel_path)

# determine path and filename
map_filename = fname.replace("_m_", "_map_")
map_filename = map_filename.replace(".fits", ".h5")
map_rel = os.path.join(sub_dir, map_filename)
# check that directory exists
if not os.path.exists(os.path.join(map_data_dir, sub_dir)):
    os.makedirs(os.path.join(map_data_dir, sub_dir), mode=0o755)

# load to LosMagneto object
full_path = os.path.join(raw_data_dir, rel_path)
hmi_im = psi_dtypes.read_hmi720s(full_path, make_map=False, solar_north_up=False)

# interpolate to map
hmi_map = hmi_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                               nprocs=nprocs, tpp=tpp, p_pool=p_pool, no_data_val=-65500.,
                               y_cor=False, helio_proj=True)

# convert interpolated map values to Br
data_index = hmi_map.data > hmi_map.no_data_val
hmi_map.data[data_index] = hmi_map.data[data_index] / hmi_map.mu[data_index]
# save Hi-Res map to file
hmi_map.y = y_axis.astype(DTypes.MAP_AXES)
hmi_map = hipft_prep.set_assim_wghts(hmi_map, assim_method="mu4_upton")
hi_res_fname = map_rel.replace("hmi_", "hmi_hires_")
hmi_map.write_to_file(map_data_dir, map_type='magneto', filename=hi_res_fname)
hmi_map.y = sin_lat.astype(DTypes.MAP_AXES)

# down-sample by integration
reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                          periodic_x=True, y_units='sinlat', uniform_poles=True,
                                          uniform_no_data=True)

# assign map y-axis back to phi
reduced_map.y = reduced_y
# set assimilation weights
reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

# write to hipft file
reduced_map.write_to_file(map_data_dir, map_type='magneto', filename=map_rel)

if nprocs > 1:
    p_pool.close()


