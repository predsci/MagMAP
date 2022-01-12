
import numpy as np
import time
import pandas as pd
import sunpy

from oftpy.utilities.datatypes import datatypes
from oftpy.maps.util import map_manip

# designate file name and location for output table
# out_path = "/Users/turtle/GitReps/OFT/oftpy/maps/devel/oftpy-interp_timings.csv"
out_path = "/Users/turtle/Dropbox/MyOFT/interp_analysis/oftpy-interp_timings.csv"
# load HMI fits file to sunpy map
fpath = "/Users/turtle/data/oft/raw_images/2014/04/13/hmi_m_720s_20140413T182401_.fits"

R0 = 1.0
out_table = pd.DataFrame(None, columns=["method", "y_pix", "map_pix", "time_s", "B_los",
                                        "B_los_abs", "B_r",	"B_r_abs", "im_min", "im_max",
                                        "im_mean", "map_min", "map_max", "map_mean"])

raw_image = sunpy.map.Map(fpath)

# --- Use OFT interp at increasing resolutions -------------------
stat_mu_thresh = 0.2

for map_nycoord in [125, 250, 500, 1000, 2000, 4000, 8000, 12000]:
    print("Starting calc for map_nycoord =", map_nycoord)
    interp_start_time = time.time()

    # the function to open the file includes the original rotation, so
    # it is included in the timing
    hmi_image = datatypes.read_hmi720s(fpath)

    map_nxcoord = 2*map_nycoord

    y_range = [-np.pi/2, np.pi/2]
    x_range = [0, 2 * np.pi]

    x_axis = np.linspace(x_range[0], x_range[1], map_nxcoord)
    y_axis = np.linspace(y_range[0], y_range[1], map_nycoord)
    # interp expects sin(lat)
    sin_lat = np.sin(y_axis)

    hmi_image.get_coordinates(R0=R0)
    hmi_image.data[np.isnan(hmi_image.data)] = 0.
    hmi_map = hmi_image.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat)
    interp_end_time = time.time()
    total_seconds = interp_end_time-interp_start_time

    print("Interp time = ", interp_end_time-interp_start_time, " seconds.")

    # Grab some comparable stats
    # ignore the limb for these stats
    stat_mu_index = hmi_image.mu >= stat_mu_thresh
    im_min = np.nanmin(raw_image.data[stat_mu_index])
    im_max = np.nanmax(raw_image.data[stat_mu_index])
    im_mean = np.nanmean(raw_image.data[stat_mu_index])

    map_data = hmi_map.data.copy()
    map_data[map_data < (hmi_map.no_data_val + 1.)] = np.nan
    stat_mu_index = hmi_map.mu >= stat_mu_thresh
    map_min = np.nanmin(map_data[stat_mu_index])
    map_max = np.nanmax(map_data[stat_mu_index])
    map_mean = np.nanmean(map_data[stat_mu_index])

    # Estimate B_los
    theta_y = np.arcsin(hmi_map.y) + np.pi/2
    theta_y[0] = 0.
    theta_y[-1] = np.pi
    phi_x = hmi_map.x
    # generate a mesh (for the accompanying area attribute)
    map_mesh = map_manip.MapMesh(phi_x, theta_y)
    # Calc B_los
    map_data[np.isnan(map_data)] = 0.
    B_data = map_mesh.da.T * map_data
    B_los = np.sum(B_data)
    B_los_abs = np.sum(np.abs(B_data))
    # Calc B_r
    data_index = (hmi_map.data > (hmi_map.no_data_val+1.)) & (hmi_map.mu > 0.05)
    B_r_data = hmi_map.data[data_index]/hmi_map.mu[data_index] * map_mesh.da.T[data_index]
    B_r_sum = np.sum(B_r_data)
    Br_abs_sum = np.sum(np.abs(B_r_data))

    # add results to output table
    out_table = out_table.append(pd.DataFrame(
        dict(method=["OFTpy_interp", ], y_pix=map_nycoord, map_pix=map_nycoord*map_nxcoord,
             time_s=total_seconds, B_los=B_los, B_los_abs=B_los_abs,
             B_r=B_r_sum, B_r_abs=Br_abs_sum, im_min=im_min, im_max=im_max,
             im_mean=im_mean, map_min=map_min, map_max=map_max, map_mean=map_mean)
    ), ignore_index=True)

# write table to csv
out_table.to_csv(out_path, index=False)
