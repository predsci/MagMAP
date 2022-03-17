
import numpy as np
import time
import pandas as pd

import sunpy
import astropy.units as u
import astropy.io
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp, reproject_adaptive, reproject_exact

from oftpy.utilities.datatypes import datatypes
from oftpy.maps.util import map_manip
from oftpy.utilities import coord_manip

# CEA fits header example
# cea_path = "/Users/turtle/Dropbox/MyOFT/fits_ex/1904-66_CEA.fits"
# map_cea = sunpy.map.Map(cea_path)
# hdulist = astropy.io.fits.open(cea_path)
# cea_header = hdulist[0].header

R0 = 1.0

y_resolutions = [125, 250, 500, 1000, 2000, 4000, 8000, 12000]
# y_resolutions = [125, 250, 500, 1000, 2000]

# load HMI fits file to sunpy map
fpath = "/Users/turtle/data/oft/raw_images/2014/04/13/hmi_m_720s_20140413T182401_.fits"
map_raw = sunpy.map.Map(fpath)
# also load into PSI format
hmi_image = datatypes.read_hmi720s(fpath)
hmi_image.get_coordinates(R0=R0)

# designate file name and location for output table
# out_path = "/Users/turtle/GitReps/OFT/oftpy/maps/devel/reproject_timings.csv"
out_path = "/Users/turtle/Dropbox/MyOFT/interp_analysis/reproject_timings.csv"

out_table = pd.DataFrame(None, columns=["method", "y_pix", "map_pix", "time_s", "B_los",
                                        "B_los_abs", "B_r",	"B_r_abs", "im_min", "im_max",
                                        "im_mean", "map_min", "map_max", "map_mean"])


def get_stats(hmi_image, raw_image, map_data, map_theta, map_phi,
              stat_mu_thresh, map_mu):
    stat_mu_index = hmi_image.mu >= stat_mu_thresh
    im_min = np.nanmin(raw_image.data[stat_mu_index])
    im_max = np.nanmax(raw_image.data[stat_mu_index])
    im_mean = np.nanmean(raw_image.data[stat_mu_index])

    # map_data[map_data < (hmi_map.no_data_val + 1.)] = np.nan
    stat_mu_index = map_mu >= stat_mu_thresh
    map_min = np.nanmin(map_data[stat_mu_index])
    map_max = np.nanmax(map_data[stat_mu_index])
    map_mean = np.nanmean(map_data[stat_mu_index])

    # Estimate B_los
    theta_y = map_theta
    theta_y[0] = 0.
    theta_y[-1] = np.pi
    phi_x = map_phi
    # generate a mesh (for the accompanying area)
    map_mesh = map_manip.MapMesh(phi_x, theta_y)
    # Calc B_los
    data_index = (~np.isnan(map_data)) & (map_mu > 0.05)
    map_data[np.isnan(map_data)] = 0.
    B_data = map_mesh.da.T*map_data
    B_los = np.sum(B_data)
    B_los_abs = np.sum(np.abs(B_data))
    # Calc B_r
    B_r_data = map_data[data_index]/map_mu[data_index]*map_mesh.da.T[data_index]
    B_r_sum = np.sum(B_r_data)
    Br_abs_sum = np.sum(np.abs(B_r_data))

    stats_out = dict(B_los=B_los, B_los_abs=B_los_abs, B_r=B_r_sum,
                     B_r_abs=Br_abs_sum, im_min=im_min, im_max=im_max,
                     im_mean=im_mean, map_min=map_min, map_max=map_max,
                     map_mean=map_mean)

    return stats_out


# --- Use reproject interp at increasing resolutions -------------------
stat_mu_thresh = 0.2

for map_nycoord in y_resolutions:
    print("Starting calc for map_nycoord =", map_nycoord)
    map_nxcoord = 2*map_nycoord

    # create Carrington WCS object
    shape_out = [map_nycoord, map_nxcoord]
    header = sunpy.map.make_fitswcs_header(np.empty(shape_out),
                                           # SkyCoord(sun.L0(map_raw.date), 0*u.deg,
                                           SkyCoord(0*u.deg,
                                                    0*u.deg,
                                                    frame="heliographic_carrington",
                                                    obstime=map_raw.date,
                                                    observer='earth'),
                                           scale=[180/map_nycoord, 360/map_nxcoord]*u.deg/u.pix,
                                           projection_code="CAR")

    # header_sinlat = sunpy.map.make_fitswcs_header(
    #     np.empty(shape_out),
    #     SkyCoord(0*u.deg, 0*u.deg, frame="heliographic_carrington",
    #              obstime=map_raw.date, observer='earth'),
    #     scale=[2*u.dimensionless_unscaled(shape_out[0]*u.pix), 0.25*u.deg/u.pix],
    #     projection_code="CEA")

    out_wcs = WCS(header)

    map_theta = np.linspace(start=0, stop=np.pi, num=map_nycoord)
    map_phi = np.linspace(start=-np.pi, stop=np.pi, num=map_nxcoord)
    sin_lat = np.sin(map_theta - np.pi/2)
    # convert 1D map axis to full list of coordinates
    mat_x, mat_y = np.meshgrid(map_phi, sin_lat)
    # convert matrix of coords to vector of coords (explicitly select row-major vectorizing)
    map_x_vec = mat_x.flatten(order="C")
    map_y_vec = mat_y.flatten(order="C")
    # clear some uneccessary vars
    mat_x = None
    mat_y = None

    # image_x, image_y, image_z, image_theta, image_phi = \
    _, _, _, image_theta, _ = \
        coord_manip.map_grid_to_image(map_x_vec, map_y_vec, R0=R0,
                                      obsv_lon=hmi_image.info['cr_lon'],
                                      obsv_lat=hmi_image.info['cr_lat'])

    # clear some uneccessary vars
    mat_x_vec = None
    mat_y_vec = None

    mu_vec = np.cos(image_theta)
    map_mu = mu_vec.reshape((map_nycoord, map_nxcoord), order="C")

    # --- reproject_interp
    interp_start_time = time.time()
    # include file-load in the timing for direct comparison
    raw_image = sunpy.map.Map(fpath)
    input_data = (raw_image.data, raw_image.wcs)

    interp_data, interp_footprint = reproject_interp(input_data, out_wcs, shape_out)

    interp_end_time = time.time()
    interp_time = interp_end_time - interp_start_time

    print("Reproject-interp time = ", interp_time, " seconds.")

    stats_dict = get_stats(hmi_image=hmi_image, raw_image=map_raw, map_data=interp_data,
                           map_theta=map_theta, map_phi=map_phi,
                           stat_mu_thresh=stat_mu_thresh, map_mu=map_mu)
    # clear large arrays
    interp_data = None
    interp_footprint = None
    # add results to output table
    out_table = out_table.append(pd.DataFrame(
        dict(method=["reproject_interp", ], y_pix=map_nycoord, map_pix=map_nycoord*map_nxcoord,
             time_s=interp_time, B_los=stats_dict["B_los"], B_los_abs=stats_dict["B_los_abs"],
             B_r=stats_dict["B_r"], B_r_abs=stats_dict["B_r_abs"], im_min=stats_dict["im_min"],
             im_max=stats_dict["im_max"], im_mean=stats_dict["im_mean"], map_min=stats_dict["map_min"],
             map_max=stats_dict["map_max"], map_mean=stats_dict["map_mean"])
        ), ignore_index=True)

    # --- reproject_adaptive
    interp_start_time = time.time()
    # the function to open the file includes the original rotation, so
    # it is included in the timing
    raw_image = sunpy.map.Map(fpath)
    input_data = (raw_image.data, raw_image.wcs)

    adaptive_data, adaptive_footprint = reproject_adaptive(input_data, out_wcs, shape_out)

    interp_end_time = time.time()
    interp_time = interp_end_time - interp_start_time

    print("Reproject-adaptive time = ", interp_time, " seconds.")

    stats_dict = get_stats(hmi_image=hmi_image, raw_image=map_raw, map_data=adaptive_data,
                           map_theta=map_theta, map_phi=map_phi,
                           stat_mu_thresh=stat_mu_thresh, map_mu=map_mu)
    # clear large arrays
    adaptive_data = None
    adaptive_footprint = None

    # add results to output table
    out_table = out_table.append(pd.DataFrame(
        dict(method=["reproject_adaptive", ], y_pix=map_nycoord, map_pix=map_nycoord*map_nxcoord,
             time_s=interp_time, B_los=stats_dict["B_los"], B_los_abs=stats_dict["B_los_abs"],
             B_r=stats_dict["B_r"], B_r_abs=stats_dict["B_r_abs"], im_min=stats_dict["im_min"],
             im_max=stats_dict["im_max"], im_mean=stats_dict["im_mean"], map_min=stats_dict["map_min"],
             map_max=stats_dict["map_max"], map_mean=stats_dict["map_mean"])
    ), ignore_index=True)

    # --- reproject_exact
    interp_start_time = time.time()
    # the function to open the file includes the original rotation, so
    # it is included in the timing
    raw_image = sunpy.map.Map(fpath)
    input_data = (raw_image.data, raw_image.wcs)

    exact_data, exact_footprint = reproject_exact(input_data, out_wcs, shape_out)

    interp_end_time = time.time()
    interp_time = interp_end_time - interp_start_time

    print("Reproject-exact time = ", interp_time, " seconds.")

    stats_dict = get_stats(hmi_image=hmi_image, raw_image=map_raw, map_data=exact_data,
                           map_theta=map_theta, map_phi=map_phi,
                           stat_mu_thresh=stat_mu_thresh, map_mu=map_mu)
    # clear large arrays
    exact_data = None
    exact_footprint = None
    # add results to output table
    out_table = out_table.append(pd.DataFrame(
        dict(method=["reproject_exact", ], y_pix=map_nycoord, map_pix=map_nycoord*map_nxcoord,
             time_s=interp_time, B_los=stats_dict["B_los"], B_los_abs=stats_dict["B_los_abs"],
             B_r=stats_dict["B_r"], B_r_abs=stats_dict["B_r_abs"], im_min=stats_dict["im_min"],
             im_max=stats_dict["im_max"], im_mean=stats_dict["im_mean"], map_min=stats_dict["map_min"],
             map_max=stats_dict["map_max"], map_mean=stats_dict["map_mean"])
    ), ignore_index=True)

    # write table to csv (at the end of each loop, in case of failure)
    out_table.to_csv(out_path, index=False)





