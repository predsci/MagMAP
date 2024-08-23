
import pandas as pd
import datetime
import os
import re

import oftpy.utilities.file_io.io_helpers as io_helpers
import oftpy.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# data-file dirs
# map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map"
map_data_dir = "/Users/turtle/data/oft/hmi_Mrmap_update"
map_type = "HMI_Mrmap"
# map_data_dir = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw"
# map_type = "HMI_Mrmap"

# movie path and frames directory
# movie_dir = "/Users/turtle/data/oft/hmi_movies/HMI_Mrmap_2012-01-15-to-2012-01-16.mov"
# movie_dir = "/Users/turtle/data/oft/hmi_movies/HMI_PSI_2012-01-15-to-2012-01-16.mov"
movie_dir = "/Users/turtle/data/oft/hmi_movies/HMI_Mrmap-update-900x900_2012-01-15-to-2012-01-16.mov"
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

f = []
for (dirpath, dirnames, filenames) in os.walk(map_data_dir):
    f.extend(filenames)
    break
available_maps = pd.DataFrame(columns=['date', 'rel_path'])
available_maps['date'] = pd.to_datetime(available_maps.date, utc=True)
available_maps.rel_path = f
available_maps = available_maps.loc[available_maps.rel_path != ".DS_Store", ]

# extract datetime from filenames
for index, row in available_maps.iterrows():
    # extract date (regex '\d{8}T\d{6}')
    date_str = re.search("\d{8}_\d{6}", row.rel_path).group()
    available_maps.loc[index, 'date'] = pd.to_datetime(date_str, format='%Y%m%d_%H%M%S', utc=True)

# ensure that maps are sorted
available_maps = available_maps.sort_values(by=['date'])

movie_df = available_maps.copy()
# change column names to 'date_mean' and 'fname'
movie_df.rename(columns=dict(date='date_mean', rel_path='fname'), inplace=True)
movie_df.reset_index(drop=True, inplace=True)

# for testing purposes, take only the first 24 maps
# movie_df = movie_df.iloc[0:24, ]

print("Converting HDF5 maps to a movie.")

psi_plt.map_movie(map_info=movie_df, png_dir=frames_dir, movie_path=movie_dir, map_dir=map_data_dir,
                  int_range=int_range, fps=fps, dpi=None, map_type=map_type, standard_grid=False)





