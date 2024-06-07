
import numpy as np

import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.utilities.plotting.psi_plotting as psi_plot
from magmap.utilities.coord_manip import *


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



hmi_path = "/Users/turtle/Dropbox/MyOFT/download_test/hmi_raw/2021/01/01/hmi_m_720s_20210101T235952_.fits"

hmi_im = psi_dtypes.read_hmi720s(hmi_path, make_map=False, solar_north_up=False)
# convert LOS B-field to B_r (with slightly larger radius to pad image for interpolation)
hmi_im.get_coordinates(R0=R0)
# convert to Br
hmi_im.add_Br(mu_thresh=0.1)

# change coordinates for interpolation
map_x = x_axis
map_y = sin_lat
interp_field = "Br"
no_data_val = -65500.
image_in = hmi_im

map_nxcoord = len(map_x)
map_nycoord = len(map_y)

# initialize grid to receive interpolation with values of NaN
interp_result = np.full((map_nycoord, map_nxcoord), no_data_val, dtype='<f4')

# convert 1D map axis to full list of coordinates
mat_x, mat_y = np.meshgrid(map_x, map_y)
# convert matrix of coords to vector of coords (explicitly select row-major vectorizing)
map_x_vec = mat_x.flatten(order="C")
map_y_vec = mat_y.flatten(order="C")
interp_result_vec = interp_result.flatten(order="C")
# determine if image is solar-north-up, or needs an additional rotation
if hasattr(image_in, "sunpy_meta"):
    if "crota2" in image_in.sunpy_meta.keys():
        image_crota2 = image_in.sunpy_meta['crota2']
    else:
        image_crota2 = 0.
else:
    image_crota2 = 0.
# convert map grid variables to image space
image_x, image_y, image_z, image_theta, image_phi = map_grid_to_image(map_x_vec, map_y_vec, R0=R0,
                                                                      obsv_lon=image_in.info['cr_lon'],
                                                                      obsv_lat=image_in.info['cr_lat'],
                                                                      image_crota2=image_crota2)
# only interpolate points on the front half of the sphere
interp_index = image_z > 0

# interpolate the specified attribute
im_data = getattr(image_in, interp_field)
interp_vec = interpolate2D_regular2irregular(image_in.x, image_in.y, im_data, image_x[interp_index],
                                             image_y[interp_index])
interp_result_vec[interp_index] = interp_vec
# reformat result to matrix form
interp_result = interp_result_vec.reshape((map_nycoord, map_nxcoord), order="C")

mu_vec = np.cos(image_theta)
mu_mat = mu_vec.reshape((map_nycoord, map_nxcoord), order="C")

out_obj = psi_dt.InterpResult(interp_result, map_x, map_y, mu_mat=mu_mat)

psi_plot.PlotDiskImage(disk_image=image_in, nfig=1, plot_attr="Br", title="Estimated B_{r} Disk")
