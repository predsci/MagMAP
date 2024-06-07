"""
Convert all maps to png in the two specified directories
"""

import os
import numpy as np

import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.utilities.file_io.io_helpers as io_helpers
import magmap.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# data-file dirs
psi_map_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800"
psi_map_cor_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_cor"
hmi_Mrmap_dir = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw"
helio_proj_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj"
helio_proj_post_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_postBr"
helio_proj_noBr_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_noBr"
helio_proj_sqrtBr_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_sqrtBr"
helio_proj_shift_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_shift"

# output dirs
png_dir = "/Users/turtle/Dropbox/MyOFT/compare-to_Mrmap/pngs_cor2"

# ----- End Inputs -------------------------

# read available magnetograms filesystem
print("\nReading filesystem from dir: " + psi_map_dir + "\n")
available_psi_maps = io_helpers.read_db_dir(psi_map_dir)

print("\nReading filesystem from dir: " + psi_map_cor_dir + "\n")
available_psi_cor_maps = io_helpers.read_db_dir(psi_map_cor_dir)

print("\nReading filesystem from dir: " + hmi_Mrmap_dir + "\n")
available_Mrmap_fits = io_helpers.read_db_dir(hmi_Mrmap_dir)

print("\nReading filesystem from dir: " + helio_proj_dir + "\n")
available_helio_proj = io_helpers.read_db_dir(helio_proj_dir)

print("\nReading filesystem from dir: " + helio_proj_post_dir + "\n")
available_helio_proj_post = io_helpers.read_db_dir(helio_proj_post_dir)

print("\nReading filesystem from dir: " + helio_proj_noBr_dir + "\n")
available_helio_proj_noBr = io_helpers.read_db_dir(helio_proj_noBr_dir)

print("\nReading filesystem from dir: " + helio_proj_sqrtBr_dir + "\n")
available_helio_proj_sqrtBr = io_helpers.read_db_dir(helio_proj_sqrtBr_dir)

print("\nReading filesystem from dir: " + helio_proj_shift_dir + "\n")
available_helio_proj_shift = io_helpers.read_db_dir(helio_proj_shift_dir)

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

# convert corrected psi_maps to png
for index, row in available_psi_cor_maps.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(psi_map_cor_dir, rel_path)
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
    png_filename = h5_date + "_psi-map_cor.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp_cor"
    # generate the png
    psi_plt.map_movie_frame(psi_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

print("\n")

# convert helio-projected psi_maps to png
for index, row in available_helio_proj.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(helio_proj_dir, rel_path)
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
    png_filename = h5_date + "_psi-helio-proj.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp_helio-proj"
    # generate the png
    psi_plt.map_movie_frame(psi_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

print("\n")

# convert helio-projected psi_maps with post-Br conversion to png
for index, row in available_helio_proj_post.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(helio_proj_post_dir, rel_path)
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
    png_filename = h5_date + "_psi-helio-proj_postBr.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp_helio-proj_postBr"
    # generate the png
    psi_plt.map_movie_frame(psi_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

print("\n")

# convert helio-projected psi_maps with no Br conversion to png
for index, row in available_helio_proj_noBr.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(helio_proj_noBr_dir, rel_path)
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
    png_filename = h5_date + "_psi-helio-proj_noBr.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp_helio-proj_noBr"
    # generate the png
    psi_plt.map_movie_frame(psi_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

print("\n")

# convert helio-projected psi_maps with no Br conversion to png
for index, row in available_helio_proj_sqrtBr.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(helio_proj_sqrtBr_dir, rel_path)
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
    png_filename = h5_date + "_psi-helio-proj_sqrtBr.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp_helio-proj_sqrtBr"
    # generate the png
    psi_plt.map_movie_frame(psi_map, int_range=[-25, 25], save_path=png_path, title=title, y_units="theta_elev",
                            no_data=True, dpi=dpi)

print("\n")

# convert helio-projected psi_maps with Yang-Br and lon-shift to png
for index, row in available_helio_proj_shift.iterrows():
    print("Converting file:", row.rel_path, "to PNG.")
    rel_path = row.rel_path
    map_path = os.path.join(helio_proj_shift_dir, rel_path)
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
    png_filename = h5_date + "_psi-helio-proj.png"
    png_path = os.path.join(png_dir, png_filename)
    # generate a title
    title = h5_date + " HMI_M-interp_helio-proj"
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

