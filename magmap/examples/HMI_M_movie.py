
import pandas as pd
import datetime
import os

import magmap.utilities.file_io.io_helpers as io_helpers
import magmap.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# data-file dirs
# map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_postBr"
map_type = "HMI_hipft"
# map_data_dir = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw"
# map_type = "HMI_Mrmap"

# movie path and frames directory
# movie_dir = "/Users/turtle/data/oft/hmi_movies/HMI_Mrmap_2012-01-15-to-2012-01-16.mov"
# movie_dir = "/Users/turtle/data/oft/hmi_movies/HMI_PSI_2012-01-15-to-2012-01-16.mov"
movie_dir = "/Users/turtle/data/oft/hmi_movies/HMI_PSI-900x900_2012-01-15-to-2012-01-16.mov"
frames_dir = "/Users/turtle/data/oft/hmi_movies/frames"

# select all maps between these dates
min_datetime_thresh = datetime.datetime(2012, 1, 15, 0, 0, 0)
max_datetime_thresh = datetime.datetime(2012, 1, 16, 0, 0, 0)

# movie frame parameters
int_range = [-25, 25]       # colorscale range
fps = 6

# ---- Main -------------------------------
# get a list of raw filenames
print("Reading map directory...")
available_maps = io_helpers.read_db_dir(map_data_dir)
# ensure that maps are sorted
available_maps = available_maps.sort_values(by=['date'])

# select only files/maps between thresh datetimes
keep_ind = (available_maps.date >= min_datetime_thresh) & \
           (available_maps.date <= max_datetime_thresh)
movie_df = available_maps.loc[keep_ind, :].copy()
# change column names to 'date_mean' and 'fname'
movie_df.rename(columns=dict(date='date_mean', rel_path='fname'), inplace=True)
movie_df.reset_index(drop=True, inplace=True)

# for testing purposes, take only the first 24 maps
# movie_df = movie_df.iloc[0:24, ]

print("Converting HDF5 maps to a movie.")

psi_plt.map_movie(map_info=movie_df, png_dir=frames_dir, movie_path=movie_dir, map_dir=map_data_dir,
                  int_range=int_range, fps=fps, dpi=None, map_type=map_type)





