
import os
import numpy as np
import sunpy

from magmap.utilities.datatypes import datatypes

# designate output map filename/path
out_path = "/Users/turtle/Dropbox/MyOFT/example_map"
out_filename = "hmi_m_720s_20140413T182401.h5"
# load HMI fits file to sunpy map
fpath = "/Users/turtle/data/oft/raw_images/2014/04/13/hmi_m_720s_20140413T182401_.fits"

R0 = 1.0
solar_north_up = False

map_nxcoord = 1024
map_nycoord = 512

# -----------------------
# query and download files

# open the file (without rotating to solar_north_up)
hmi_image = datatypes.read_hmi720s(fpath, solar_north_up=solar_north_up)
hmi_image.get_coordinates(R0=R0)
hmi_image.data[np.isnan(hmi_image.data)] = 0.

# setup map grid
y_range = [-np.pi/2, np.pi/2]
x_range = [0, 2 * np.pi]
x_axis = np.linspace(x_range[0], x_range[1], map_nxcoord)
y_axis = np.linspace(y_range[0], y_range[1], map_nycoord)
# interp expects sin(lat)
sin_lat = np.sin(y_axis)

# interp to map
hmi_map = hmi_image.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat)

# assign map y-axis back to phi
hmi_map.y = y_axis

# write to hipft file
hmi_map.write_to_file(out_path, map_type='magneto', filename=out_filename)


test = datatypes.read_hipft_map(os.path.join(out_path, out_filename))

raw_image = sunpy.map.Map(fpath)
raw_image.peek()
