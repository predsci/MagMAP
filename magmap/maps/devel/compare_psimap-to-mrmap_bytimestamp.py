"""
Given an image relative path, generate maps by three different methods:
  - JSOC method to produce HMI_Mrmap
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
# image relative path and timestamp string
im_fits_dir = "/Volumes/extdata3/oft/raw_data/hmi_m720s"
im_rel_path = "2012/01/15/hmi_m_720s_20120115T015952.fits"
im_timestamp = "20120115T015952"

# Mrmap file dirs
hmi_Mrmap_dir = "/Volumes/extdata3/oft/raw_data/hmi_mrmap_720s"
Mrmap_rel_path = "2012/01/15/hmi_Mrmap_latlon_720s_20120115T015952.fits"

# input full paths
hmi_m_fits = os.path.join(im_fits_dir, im_rel_path)
hmi_mram_fits = os.path.join(hmi_Mrmap_dir, Mrmap_rel_path)

# output dirs
png_dir = "/Users/turtle/Dropbox/MyOFT/compare-to_Mrmap/pngs_shift"
write_maps = True
map_data_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_shift"
psi_shifted_dir = "/Volumes/terminus_ext/HMI_M720s/hmi_map_900x1800_helio-proj_shift"

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

# load the image
hmi_m_im = psi_dtypes.read_hmi720s(hmi_m_fits)
# load corresponding HMI_Mrmap
mr_map = psi_dtypes.read_hmi_Mrmap_latlon_720s(hmi_mram_fits)

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
reduced_x = np.linspace(0.1, 359.9, reduced_nxcoord) * np.pi/180
reduced_y = np.linspace(-89.9, 89.9, reduced_nycoord) * np.pi/180
# interp expects sin(lat)
reduced_sin_lat = np.sin(reduced_y)

# generate filename
fname = os.path.basename(im_rel_path)
sub_dir = os.path.dirname(im_rel_path)

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

# record Mrmap pixels for chopping shifted images
# Mrmap_data_index = mr_map.data > mr_map.no_data_val
Mrmap_data_index = (mr_map.data > mr_map.no_data_val)[1:-1, 1:-1]
Mrmap_data_vals = mr_map.data[1:-1, 1:-1][Mrmap_data_index]

# loop through a series of left-right pixel shifts
n01_shift_min = -20
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
plt.title(im_timestamp)
plt.savefig(png_dir + "/median_diff_" + im_timestamp + ".png", dpi=600)
plt.close()


plt.plot(shift_deg, diff_avg)
# plt.plot(shift_deg, diff_rms)
plt.xlabel("Longitudinal shift (deg)")
plt.ylabel("Mean absolute difference (Gs)")
# plt.ylim([0., 6])
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.title(im_timestamp)
plt.savefig(png_dir + "/mean_diff_" + im_timestamp + ".png", dpi=600)
plt.close()


plt.plot(shift_deg, diff_rms)
plt.xlabel("Longitudinal shift (deg)")
plt.ylabel("RMS difference (Gs)")
# plt.ylim([0., 6])
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.title(im_timestamp)
plt.savefig(png_dir + "/rms_diff_" + im_timestamp + ".png", dpi=600)
plt.close()

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

# auto-choose best shift?
shifts = 0.1*(np.array(range(n01_shift_min, n01_shift_max))/10)
best_shift = np.round(100*shifts[np.argmin(diff_median)])
shift_str = f'{abs(best_shift):02.0f}'
if best_shift < 0.:
    sign_str = "m"
elif best_shift > 0.:
    sign_str = "p"
else:
    sign_str = ""
shifted_rel = map_rel.replace(".", "_shifted" + sign_str + shift_str + ".")

# set map path and open
best_shift_path = os.path.join(map_data_dir, shifted_rel)
best_shift_map = psi_dtypes.read_hipft_map(best_shift_path)

# best-shift with PSI Br
best_shift_path_psi = os.path.join(psi_shifted_dir, shifted_rel)
best_shift_psi = psi_dtypes.read_hipft_map(best_shift_path_psi)

# psi_data = np.copy(best_shift_map.data[1:-1, 1:-1])
psi_data = np.copy(best_shift_map.data)
jitter_data = np.copy(mr_map.data[1:-1, 1:-1])

# convert jitter-mu back to correct mu
jitter_mu = np.copy(best_shift_map.mu)
jitter_im_radius = np.sqrt(1 - jitter_mu**2)
radius_per_half_image = r_sun_obs/cdelt/im_pix_width
im_radius = jitter_im_radius/radius_per_half_image
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
jitter_mean = np.full([len(mu_edge_vec)-1, ], 0.)
psi_abs_mean = np.full([len(mu_edge_vec)-1, ], 0.)
jitter_abs_mean = np.full([len(mu_edge_vec)-1, ], 0.)
psiBr_mean = np.full([len(mu_edge_vec)-1, ], 0.)
psiBr_abs_mean = np.full([len(mu_edge_vec)-1, ], 0.)

# PSI maps tend to have slightly less pixels at edge because of edge-preserving downsample
psi_index = psi_data > best_shift_map.no_data_val
jitter_index = jitter_data > mr_map.no_data_val
psiBr_index = best_shift_psi.data != best_shift_psi.no_data_val
data_index = psi_index & psiBr_index & jitter_index
for ii in range(len(solar_radii_vec)-1):
    # determine which pixels land in this bin
    upper_thresh_mu = mu_edge_vec[ii]
    lower_thresh_mu = mu_edge_vec[ii+1]
    mu_index = (cor_mu <= upper_thresh_mu) & (cor_mu > lower_thresh_mu)
    best_shift_index = mu_index & data_index
    # calculate a mean value for each bin
    psi_mean[ii] = np.nanmean(psi_data[best_shift_index])
    jitter_mean[ii] = np.nanmean(jitter_data[best_shift_index])
    psi_abs_mean[ii] = np.nanmean(np.abs(psi_data[best_shift_index]))
    jitter_abs_mean[ii] = np.nanmean(np.abs(jitter_data[best_shift_index]))
    # also record for PSI-Br
    psiBr_mean[ii] = np.nanmean(best_shift_psi.data[best_shift_index])
    psiBr_abs_mean[ii] = np.nanmean(np.abs(best_shift_psi.data[best_shift_index]))

center_limb_ang = np.arcsin(solar_radii_centers)*180/np.pi

# first plot PSI-Br vs jitter Br
plt.figure(1)
plt.plot(center_limb_ang, psiBr_abs_mean)
plt.plot(center_limb_ang, jitter_abs_mean)
plt.yscale('log')
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("Mean Abs magnitude (Gs)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.legend(['PSI', 'HMI_Mrmap'])

plt.savefig(png_dir + "/Br-PSI-jitter_v_radius_" + im_timestamp + ".png", dpi=600)
plt.close(1)


# first plot PSI-Br div jitter Br
plt.figure(1)
plt.plot(deg_vec[:-11], los2br[:-11]/los2br_cor[:-11])
plt.plot(center_limb_ang[:-1], psiBr_abs_mean[:-1]/jitter_abs_mean[:-1])
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("Mean Abs magnitude (PSI/HMI_Mrmap)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)
plt.legend(["Theoretic", "Actual"])

plt.savefig(png_dir + "/Br-PSI-div-jitter_v_radius_" + im_timestamp + ".png", dpi=600)
plt.close(1)


# plot the ratio between absolute means
plt.figure(1)
plt.plot(center_limb_ang, psi_abs_mean/jitter_abs_mean)
plt.xlabel("Center-to-limb angle (deg)")
plt.ylabel("Mean Abs (PSI/jitter)")
plt.grid(alpha=0.6, linestyle='dashed', lw=0.5)

plt.savefig(png_dir + "/PSI-div-jitter_v_radius_" + im_timestamp + ".png", dpi=600)
plt.close(1)
