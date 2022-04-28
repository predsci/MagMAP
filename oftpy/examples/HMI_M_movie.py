
import pandas as pd
import datetime
import os

import oftpy.utilities.file_io.io_helpers as io_helpers
import oftpy.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# data-file dirs
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"

# movie path and frames directory
movie_dir = "/Users/turtle/data/oft/hmi_movies/2011-02-20-to-2011-04-01.mov"
frames_dir = "/Users/turtle/data/oft/hmi_movies/frames"

# select all maps between these dates
min_datetime_thresh = datetime.datetime(2011, 2, 20, 0, 0, 0)
max_datetime_thresh = datetime.datetime(2011, 4, 1, 0, 0, 0)

# movie frame parameters
int_range = [-25, 25]       # colorscale range


# ---- Main -------------------------------
# get a list of raw filenames
print("Reading map directory...")
available_maps = io_helpers.read_db_dir(map_data_dir)

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
                  int_range=int_range, fps=24, dpi=None, map_type="HMI_hipft")





