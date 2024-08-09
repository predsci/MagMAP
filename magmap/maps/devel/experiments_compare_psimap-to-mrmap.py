"""
Given an image, generate maps by three different methods:
  - Yang/JSOC method to produce HMI_Mrmap
  - PSI method to reproduce HMI_Mrmap
  - PSI method
Then do multiple comparisons of the LOS-to-Br conversion and longitudidnal shift
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

import magmap.utilities.datatypes.datatypes as psi_dtypes
import magmap.maps.hipft_prep as hipft_prep
import magmap.maps.util.map_manip as map_manip
# additional routines for debug/development
import magmap.utilities.file_io.io_helpers as io_helpers
import magmap.utilities.plotting.psi_plotting as psi_plt

# ---- Inputs -----------------------------
# data-file dirs
psi_map_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800"
hmi_Mrmap_dir = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw"

# output dirs
png_dir = "/Users/turtle/Dropbox/MyOFT/compare-to_Mrmap/pngs_shift"
write_maps = False
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_shift"
psi_shifted_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_shift"

# fits paths
hmi_m_fits = "/Volumes/terminus_ext/HMI_M720s/hmi_raw/2012/01/15/hmi_m_720s_20120115T005952.fits"
hmi_mram_fits = "/Volumes/terminus_ext/HMI_Mrmap_latlon_720s/hmi_raw/2012/01/15/" + \
                "hmi_Mrmap_latlon_720s_20120115T005952.fits"
# map files
psi_cur = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_postBr/2012/01/15/hmi_map_720s_20120115T005952.h5"
psi_sqrtBr = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_sqrtBr/2012/01/15/hmi_map_720s_20120115T005952.h5"

# number of processors for interpolation parallelization
nprocs = 4
tpp = 5

# High-res map grid specifications
map_nxcoord = 10240
map_nycoord = 5120
R0 = 1.

# reduced map grid specifications
reduced_nxcoord = 1800
reduced_nycoord = 900

# ----- End Inputs -------------------------

hmi_m_im = psi_dtypes.read_hmi720s(hmi_m_fits)

mr_map = psi_dtypes.read_hmi_Mrmap_latlon_720s(hmi_mram_fits)

psi_map = psi_dtypes.read_hipft_map(psi_cur)

mram_index = mr_map.data > mr_map.no_data_val
test = psi_map.data / mr_map.data

test[~mram_index] = np.nan

plt.imshow(test, vmin=-2, vmax=2)
plt.colorbar()


psi_sqrtBr_map = psi_dtypes.read_hipft_map(psi_sqrtBr)
test = psi_sqrtBr_map.data / mr_map.data
test[~mram_index] = np.nan

plt.imshow(test, vmin=-2, vmax=2)
plt.colorbar()


## --- test map shift, compare to Mrmap ----------------------

# initiate a pool of processors
p_pool = Pool(nprocs)

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

# setup reduced-map grid (specific to hmi_Mrmap_latlon_720s grid),
# with added points at poles and periodic boundary to match OFTpy
# mapping standard
reduced_x = np.linspace(0.1, 359.9, reduced_nxcoord) * np.pi/180
reduced_y = np.linspace(-89.9, 89.9, reduced_nycoord) * np.pi/180
# interp expects sin(lat)
reduced_sin_lat = np.sin(reduced_y)

# generate filename
rel_path = '2012/01/15/hmi_m_720s_20120115T005952.fits'
fname = os.path.basename(rel_path)
sub_dir = os.path.dirname(rel_path)

# determine path and filename
map_filename = fname.replace("_m_", "_map_")
map_filename = map_filename.replace(".fits", ".h5")
map_rel = os.path.join(sub_dir, map_filename)

# check that directory exists
if not os.path.exists(os.path.join(map_data_dir, sub_dir)):
    os.makedirs(os.path.join(map_data_dir, sub_dir), mode=0o755)

hmi_map = hmi_m_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                 nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=True)
hmi_map.y = y_axis

# convert interpolated map values to Br
data_index = hmi_map.data > hmi_map.no_data_val
hmi_map.data[data_index] = hmi_map.data[data_index] / hmi_map.mu[data_index]

# down-sample by integration (normal grid)
reduced_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=reduced_x, image_method=0,
                                          periodic_x=True, y_units='sinlat', uniform_poles=True,
                                          uniform_no_data=True)

# also interpolate by PSI method (to calc Br PSI method)
psi_map = hmi_m_im.interp_to_map(R0=R0, map_x=x_axis, map_y=sin_lat, interp_field="data",
                                 nprocs=nprocs, tpp=tpp, p_pool=p_pool, y_cor=False, helio_proj=True)
psi_map.y = y_axis

# convert interpolated map values to Br
data_index = psi_map.data > psi_map.no_data_val
psi_map.data[data_index] = psi_map.data[data_index] / psi_map.mu[data_index]

# check that directory exists
if not os.path.exists(os.path.join(psi_shifted_dir, sub_dir)):
    os.makedirs(os.path.join(psi_shifted_dir, sub_dir), mode=0o755)

# assign map y-axis back to phi
reduced_map.y = reduced_y
# set assimilation weights
reduced_map = hipft_prep.set_assim_wghts(reduced_map, assim_method="mu4_upton")

# write to hipft file
reduced_map.write_to_file(map_data_dir, map_type='magneto', filename=map_rel)

# record Mrmap pixels for chopping shifted images
# Mrmap_data_index = mr_map.data > mr_map.no_data_val
Mrmap_data_index = (mr_map.data > mr_map.no_data_val)[1:-1, 1:-1]
Mrmap_data_vals = mr_map.data[1:-1, 1:-1][Mrmap_data_index]

# loop through a series of left-right pixel shifts
n01_shift_min = -10
n01_shift_max = 30 + 1
nshifts = n01_shift_max - n01_shift_min
diff_avg = np.full([nshifts, ], fill_value=0.)
diff_median = np.full([nshifts, ], fill_value=0.)
diff_rms = np.full([nshifts, ], fill_value=0.)
for ii in range(n01_shift_min, n01_shift_max):
    rad_shift = 0.1*np.pi/180 * (ii/10)
    # generate a new longitudinal axis
    shifted_x = reduced_x + rad_shift
    left_pix_condense = shifted_x <= 0.
    n_left_con_pix = np.sum(left_pix_condense)
    if n_left_con_pix > 0.:
        # condense pixels at left edge
        left_margin = shifted_x[n_left_con_pix]
        shifted_x[left_pix_condense] = left_margin*(np.argwhere(left_pix_condense)[:, 0] + 0.5)/n_left_con_pix
        # reset right margin to 2pi (?)
        shifted_x[-1] = 2 * np.pi
    else:
        right_pix_condense = shifted_x >= 2*np.pi
        n_right_con_pix = np.sum(right_pix_condense)
        if n_right_con_pix > 0.:
            # condense pixels at right edge
            right_interior = shifted_x[-(n_right_con_pix+1)]
            right_margin = 2*np.pi - right_interior
            right_rel_index = (shifted_x.__len__() - 1) - np.argwhere(right_pix_condense)[:, 0]
            shifted_x[right_pix_condense] = 2*np.pi - right_margin * (right_rel_index + 0.5)/(n_right_con_pix)
            # reset left margin to 0 (?)
            shifted_x[0] = 0

    print("Calculating map shifted by", f'{ii*0.01:3.2}', "degrees longitude.")
    # down-sample by integration (shifted grid)
    # shifted_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_sin_lat, new_x=shifted_x, image_method=0,
    #                                           periodic_x=True, y_units='sinlat', uniform_poles=True,
    #                                           uniform_no_data=True)
    shifted_map = map_manip.downsamp_reg_grid(full_map=hmi_map, new_y=reduced_y, new_x=shifted_x,
                                              new_x_extents=x_range, new_y_extents=y_range,
                                              image_method=0, periodic_x=False, y_units='lat_rad',
                                              uniform_poles=False, uniform_no_data=True)
    shifted_map.data[~Mrmap_data_index] = shifted_map.no_data_val
    # assign map y-axis back to phi
    shifted_map.y = reduced_y
    if write_maps:
        # set assimilation weights
        shifted_map = hipft_prep.set_assim_wghts(shifted_map, assim_method="mu4_upton")

        if ii < 0:
            sign_str = "m"
        elif ii > 0.:
            sign_str = "p"
        else:
            sign_str = ""
        shifted_rel = map_rel.replace(".", "_shifted" + sign_str + f'{abs(ii):02}' + ".")
        # write to hipft file
        shifted_map.write_to_file(map_data_dir, map_type='magneto', filename=shifted_rel)

    # compare to Mrmap
    shifted_data_vals = shifted_map.data[Mrmap_data_index]
    map_diff = np.abs(shifted_data_vals - Mrmap_data_vals)
    vector_index = ii - n01_shift_min
    diff_avg[vector_index] = np.sum(map_diff)/Mrmap_data_vals.__len__()
    diff_median[vector_index] = np.median(map_diff)
    diff_rms[vector_index] = np.sqrt(np.sum(map_diff**2))/Mrmap_data_vals.__len__()

    # also create a shifted map using PSI-Br method
    shifted_map = map_manip.downsamp_reg_grid(full_map=psi_map, new_y=reduced_y, new_x=shifted_x,
                                              new_x_extents=x_range, new_y_extents=y_range,
                                              image_method=0, periodic_x=False, y_units='lat_rad',
                                              uniform_poles=False, uniform_no_data=True)
    shifted_map.data[~Mrmap_data_index] = shifted_map.no_data_val
    if write_maps:
        # set assimilation weights
        shifted_map = hipft_prep.set_assim_wghts(shifted_map, assim_method="mu4_upton")
        # write to hipft file
        shifted_map.write_to_file(psi_shifted_dir, map_type='magneto', filename=shifted_rel)

shift_deg = np.multiply(list(range(n01_shift_min, n01_shift_max)), 0.1)*0.1

plt.plot(shift_deg, diff_median)
# plt.plot(shift_deg, diff_avg)
# plt.plot(shift_deg, diff_rms)
plt.xlabel("Longitudinal shift (deg)")
plt.ylabel("Median absolute difference (Gs)")
plt.ylim([0., 6])
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.title("2012/01/15 00:59:52")
plt.savefig(png_dir + "/median_diff_20120115T005952.png", dpi=600)
plt.close()


plt.plot(shift_deg, diff_avg)
# plt.plot(shift_deg, diff_rms)
plt.xlabel("Longitudinal shift (deg)")
plt.ylabel("Mean absolute difference (Gs)")
# plt.ylim([0., 6])
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.title("2012/01/15 00:59:52")
plt.savefig(png_dir + "/mean_diff_20120115T005952.png", dpi=600)
plt.close()


plt.plot(shift_deg, diff_rms)
plt.xlabel("Longitudinal shift (deg)")
plt.ylabel("RMS difference (Gs)")
# plt.ylim([0., 6])
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.title("2012/01/15 00:59:52")
plt.savefig(png_dir + "/rms_diff_20120115T005952.png", dpi=600)

r_sun_obs = hmi_m_im.sunpy_meta['rsun_obs']
cdelt = hmi_m_im.sunpy_meta['cdelt1']
im_pix_width = hmi_m_im.data.shape[0]/2

theta_vec = np.linspace(start=0., stop=np.pi/2, num=100)
deg_vec = theta_vec*180/np.pi
mu_vec = np.cos(theta_vec)
los2br = 1/mu_vec

sin_vec = np.sin(theta_vec)
sin_vec_cor = (r_sun_obs/cdelt/im_pix_width) * sin_vec
los2br_cor = 1/np.sqrt(1 - sin_vec_cor**2)

# plot losbr and los2br_cor on same plot
plt.figure(1)
plt.plot(deg_vec[:-1], los2br[:-1])
plt.plot(deg_vec[:-1], los2br_cor[:-1])
plt.yscale("log")
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("LOS-to-Br multiplicative factor")
plt.legend(["PSI", "Mrmap"])
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.savefig(png_dir + "/los2br_both.png", dpi=600)
plt.close(1)

# plot losbr/los2br
plt.figure(1)
plt.plot(deg_vec[:-11], los2br[:-11]/los2br_cor[:-11])
# plt.yscale("log")
plt.title("LOS-to-Br: PSI v Mrmap")
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("PSI/Mrmap")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.savefig(png_dir + "/los2br_div.png", dpi=600)
plt.close(1)

# --- re-plot on a sinlat axis ------
theta_vec2 = np.linspace(start=-np.pi/2, stop=np.pi/2, num=200)
deg_vec2 = theta_vec2*180/np.pi
mu_vec2 = np.cos(theta_vec2)
los2br_2 = 1/mu_vec2

sin_vec2 = np.sin(theta_vec2)
sin_vec_cor2 = (r_sun_obs/cdelt/im_pix_width) * sin_vec2
los2br_cor_2 = 1/np.sqrt(1 - sin_vec_cor2**2)
# plot losbr/los2br
plt.figure(1)
plt.plot(sin_vec2[6:-6], los2br_2[6:-6]/los2br_cor_2[6:-6])
# plt.yscale("log")
plt.title("LOS-to-Br: PSI v Mrmap")
plt.xlabel("Sine Lat")
plt.ylabel("PSI/Mrmap")
plt.ylim((-2, 4))
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.savefig(png_dir + "/los2br_div_sinlat.png", dpi=600)
plt.close(1)


# plot losbr - los2br
plt.figure(1)
plt.plot(deg_vec[:-11], los2br[:-11] - los2br_cor[:-11])
# plt.yscale("log")
plt.title("LOS-to-Br: PSI v Mrmap")
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("PSI - Mrmap")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.savefig(png_dir + "/los2br_diff.png", dpi=600)
plt.close(1)


## --- compare average values by radius -----------------
# set map path and open
best_shift_path = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_shift/2012/01/15/" \
                  "hmi_map_720s_20120115T005952_shiftedp13.h5"
# best_shift_path = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_shift/2012/01/15/" \
#                   "hmi_map_720s_20120115T005952_shifted00.h5"
best_shift_map = psi_dtypes.read_hipft_map(best_shift_path)

# best-shift with PSI Br
best_shift_psi = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_shift/2012/01/15/" \
                  "hmi_map_720s_20120115T005952_shiftedp13.h5"
best_shift_psi = psi_dtypes.read_hipft_map(best_shift_psi)

psi_data = np.copy(best_shift_map.data[1:-1, 1:-1])
yang_data = np.copy(mr_map.data[1:-1, 1:-1])

# convert yang-mu back to correct mu
yang_mu = np.copy(best_shift_map.mu[1:-1, 1:-1])
yang_im_radius = np.sqrt(1 - yang_mu**2)
radius_per_half_image = r_sun_obs/cdelt/im_pix_width
im_radius = yang_im_radius/radius_per_half_image
cor_mu = np.sqrt(1 - im_radius**2)

# determine observer CR coords from fits header
CR_lat = hmi_m_im.sunpy_meta["crlt_obs"]
CR_lon = hmi_m_im.sunpy_meta["crln_obs"]
# set radius bin edges
solar_radii_vec = np.sqrt(np.linspace(start=0, stop=1, num=41))
solar_radii_centers = (solar_radii_vec[:-1] + solar_radii_vec[1:])/2
mu_edge_vec = np.sqrt(1 - solar_radii_vec**2)
mu_center_vec = (mu_edge_vec[:-1] + mu_edge_vec[1:])/2
psi_mean = np.full([len(mu_edge_vec)-1, ], 0.)
yang_mean = np.full([len(mu_edge_vec)-1, ], 0.)
psi_abs_mean = np.full([len(mu_edge_vec)-1, ], 0.)
yang_abs_mean = np.full([len(mu_edge_vec)-1, ], 0.)
psiBr_mean = np.full([len(mu_edge_vec)-1, ], 0.)
psiBr_abs_mean = np.full([len(mu_edge_vec)-1, ], 0.)

# PSI maps tend to have slightly less pixels at edge because of edge-preserving downsample
psi_index = psi_data > best_shift_map.no_data_val
yang_index = yang_data > mr_map.no_data_val
psiBr_index = best_shift_psi.data != best_shift_psi.no_data_val
data_index = psi_index & psiBr_index & yang_index
for ii in range(len(solar_radii_vec)-1):
    # determine which pixels land in this bin
    upper_thresh_mu = mu_edge_vec[ii]
    lower_thresh_mu = mu_edge_vec[ii+1]
    mu_index = (cor_mu <= upper_thresh_mu) & (cor_mu > lower_thresh_mu)
    best_shift_index = mu_index & data_index
    # calculate a mean value for each bin
    psi_mean[ii] = np.nanmean(psi_data[best_shift_index])
    yang_mean[ii] = np.nanmean(yang_data[best_shift_index])
    psi_abs_mean[ii] = np.nanmean(np.abs(psi_data[best_shift_index]))
    yang_abs_mean[ii] = np.nanmean(np.abs(yang_data[best_shift_index]))
    # also record for PSI-Br
    psiBr_mean[ii] = np.nanmean(best_shift_psi.data[best_shift_index])
    psiBr_abs_mean[ii] = np.nanmean(np.abs(best_shift_psi.data[best_shift_index]))

center_limb_ang = np.arcsin(solar_radii_centers)*180/np.pi

# first plot PSI-Br vs Yang Br
plt.figure(0)
plt.plot(center_limb_ang, psiBr_abs_mean)
plt.plot(center_limb_ang, yang_abs_mean)
plt.yscale('log')
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("Mean Abs magnitude (Gs)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.legend(['PSI', 'HMI_Mrmap'])

plt.savefig(png_dir + "/Br-PSI-Yang_v_radius.png", dpi=600)
plt.close(0)


# first plot PSI-Br div Yang Br
plt.figure(0)
plt.plot(deg_vec[:-11], los2br[:-11]/los2br_cor[:-11])
plt.plot(center_limb_ang[:-1], psiBr_abs_mean[:-1]/yang_abs_mean[:-1])
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("Mean Abs magnitude (PSI/HMI_Mrmap)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.legend(["Theoretic", "Actual"])

plt.savefig(png_dir + "/Br-PSI-div-Yang_v_radius.png", dpi=600)
plt.close(0)


# plot the ratio between absolute means
plt.figure(0)
plt.plot(center_limb_ang, psi_abs_mean/yang_abs_mean)
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("Mean Abs (PSI/Yang)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)

plt.savefig(png_dir + "/PSI-div-Yang_v_radius.png", dpi=600)
plt.close(0)


# determine PSI/MRmap ratio for each Theta value
psi_data = np.copy(best_shift_map.data)
ratio_data = np.abs(psi_data/yang_data)

ignore_index = (np.isnan(ratio_data)) | (ratio_data > 1000)
row_means = np.mean(ratio_data, axis=1, where=~ignore_index)
row_means2 = np.mean(ratio_data, axis=1)

theta_lat = mr_map.y[1:-1]
im_shifted_theta = theta_lat - CR_lat*np.pi/180
use_index = (im_shifted_theta >= -np.pi/2) & (im_shifted_theta <= np.pi/2)
im_sine_lat = np.sin(im_shifted_theta)

# plot the ratio between methods as a function of sine_lat
plt.figure(0)
plt.plot(im_sine_lat[use_index], row_means[use_index])
plt.xlabel("Sine Lat")
plt.ylabel("Mean Abs (PSI/Yang)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)

plt.savefig(png_dir + "/PSI-div-Yang_v_radius.png", dpi=600)
plt.close(0)

# determine PSI/MRmap ratio for each Theta value (theoretical)
yang_mu = np.copy(best_shift_map.mu)
yang_im_radius = np.sqrt(1 - yang_mu**2)
radius_per_half_image = r_sun_obs/cdelt/im_pix_width
im_radius = yang_im_radius/radius_per_half_image
psi_mu = np.sqrt(1 - im_radius**2)

multi_factor = yang_mu/psi_mu

psi_index = best_shift_map.data > best_shift_map.no_data_val
yang_index = mr_map.data[1:-1, 1:-1] > mr_map.no_data_val
psiBr_index = best_shift_psi.data != best_shift_psi.no_data_val
data_index = psi_index & psiBr_index & yang_index

mu_limit = 0.1
keep_index = (psi_mu > mu_limit) & data_index
row_means = np.mean(multi_factor, axis=1, where=keep_index)

theta_lat = mr_map.y[1:-1]
im_shifted_theta = theta_lat - CR_lat*np.pi/180
use_index = (im_shifted_theta >= -np.pi/2) & (im_shifted_theta <= np.pi/2)
im_sine_lat = np.sin(im_shifted_theta)

# plot the ratio between methods as a function of sine_lat
plt.figure(0)
plt.plot(im_sine_lat[use_index], row_means[use_index])
plt.ylim((-2, 4))
plt.xlabel("Sine Lat")
plt.ylabel("Mean (PSI/Yang)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)

plt.savefig(png_dir + "/PSI-div-Yang_v_sinlat_" + str(mu_limit) + ".png", dpi=600)
plt.close(0)

