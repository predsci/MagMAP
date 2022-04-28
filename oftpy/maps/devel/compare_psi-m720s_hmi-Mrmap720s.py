"""
Convert all maps to png in the two specified directories
"""

import os
import numpy as np

import oftpy.utilities.datatypes.datatypes as psi_dtypes
import oftpy.utilities.file_io.io_helpers as io_helpers
import oftpy.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# data-file dirs
psi_map_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800"
hmi_Mrmap_dir = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw"

# output dirs
png_dir = "/Users/turtle/Dropbox/MyOFT/compare-to_Mrmap/pngs"

# ----- End Inputs -------------------------

# read available magnetograms filesystem
print("\nReading filesystem from dir: " + psi_map_dir + "\n")
available_psi_maps = io_helpers.read_db_dir(psi_map_dir)

print("\nReading filesystem from dir: " + hmi_Mrmap_dir + "\n")
available_Mrmap_fits = io_helpers.read_db_dir(hmi_Mrmap_dir)


# convert psi_maps to png
for index, row in available_psi_maps.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(psi_map_dir, rel_path)
    psi_map = psi_dtypes.read_hipft_map(map_path)

    # determine 'full' dpi
    map_shape = psi_map.data.shape
    width_dpi = map_shape[1] / 5.57
    height_dpi = map_shape[0] / 2.78
    dpi = np.ceil(max(width_dpi, height_dpi))

    # generate a filename leading with date
    h5_filename = os.path.basename(rel_path)
    h5_date = h5_filename.replace("hmi_map_720s_", "")
    h5_date = h5_date.replace(".h5", "")
    png_filename = h5_date + "_psi-map.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp"
    # generate the png
    psi_plt.map_movie_frame(psi_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

print("\n")

# convert Mrmaps to png
for index, row in available_Mrmap_fits.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(hmi_Mrmap_dir, rel_path)
    Mrmap_map = psi_dtypes.read_hmi_Mrmap_latlon_720s(map_path)

    # determine 'full' dpi
    map_shape = Mrmap_map.data.shape
    width_dpi = map_shape[1] / 5.57
    height_dpi = map_shape[0] / 2.78
    dpi = np.ceil(max(width_dpi, height_dpi))

    # generate a filename leading with date
    fits_filename = os.path.basename(rel_path)
    fits_date = fits_filename.replace("hmi_Mrmap_latlon_720s_", "")
    fits_date = fits_date.replace(".fits", "")
    png_filename = fits_date + "_Mrmap.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = fits_date + " HMI_Mrmap"
    # generate the png
    psi_plt.map_movie_frame(Mrmap_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

