"""
Functions to manipulate and combine maps
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy import sparse
import time

import magmap.utilities.datatypes.datatypes as psi_d_types
from magmap.utilities.coord_manip import s2c


class MapMesh:
    """Simple class to hold mesh information for a phi/theta map.

    The purpose of this class is to facilitate map/mesh based computations such as integral
    averaging or mesh spacing calculations.

    This holds things like:
    - the exterior (main) mesh locations.
    - the interior (half) mesh locations.
    - the area factors at the main mesh points (using an approximation for the edges).
    - 2D arrays with the cartesian locations of each point (for easy averaging calculations).
    - interpolators to get the cell index given a phi or theta value
    - interpolators to get the local mesh spacing for a given phi or theta value

    The phi and theta can be non-uniformly spaced.

    This class is not as flexible/general as what is in our fortran code remesh_general
    but it is suitable for our purposes.

    We can add other interpolators / interpolation methods for maps to this later if needed.
    """

    def __init__(self, p, t):
        """
        Generate the mesh object based on 1D arrays of phi and theta locations
        """
        # define a floating point tolerance for checking if at 0, pi, or 2pi
        eps = 2e-7

        # check to see if it is periodic (ASSUME ONE REPEATED POINT)
        periodic = False
        if abs(p[0]) <= eps and abs(p[-1] - 2*np.pi) <= eps:
            periodic = True

        # get the 1D mesh properties for phi and theta:
        #   this is: half mesh locations, main mesh spacing, half mesh spacing
        ph, dp, dph = get_1d_mesh_properties(p)
        th, dt, dth = get_1d_mesh_properties(t)

        # get a special version of dp for periodic interpolation (put the whole cell width vs. half)
        # need 2 versions since for area, you want the half cell!
        dp_alt = np.copy(dp)
        if periodic:
            dp_alt[0] = dp[0] + dp[-1]
            dp_alt[-1] = dp_alt[0]

        # Get the coordinates, assuming radius = 1 (photosphere)
        # here i want to make an array that is phi, theta ordered on purpose
        t2d, p2d = np.meshgrid(t, p)
        dt2d, dp2d = np.meshgrid(dt, dp)
        r2d = np.ones_like(t2d)

        # get the cartesian coordinates (for convenience in later computations)
        x2d, y2d, z2d = s2c(r2d, t2d, p2d)

        # get the area of each pixel, but modify t2d to account for polar points (shift a quarter point)
        if abs(t[0]) <= eps:
            t2d[:, 0] = 0.5*(t[0] + th[0])
        if abs(t[-1] - np.pi) <= eps:
            t2d[:, -1] = 0.5*(t[-1] + th[-1])
        da = np.sin(t2d)*dt2d*dp2d

        # now build interpolators that will be helpful
        interp_p2index = interp1d(p, np.arange(len(p)), fill_value=(0, len(p) - 1), bounds_error=False)
        interp_t2index = interp1d(t, np.arange(len(t)), fill_value=(0, len(t) - 1), bounds_error=False)
        interp_p2dp = interp1d(p, dp_alt, fill_value=(dp_alt[0], dp_alt[-1]), bounds_error=False)
        interp_t2dt = interp1d(t, dt, fill_value=(dt[0], dt[-1]), bounds_error=False)

        # add these as attributes to the class
        self.n_p = len(p)
        self.n_t = len(t)
        self.p = p
        self.t = t
        self.dp = dp
        self.dt = dt
        self.ph = ph
        self.th = th
        self.dph = dph
        self.dth = dth
        self.da = da
        self.periodic = periodic
        self.interp_p2index = interp_p2index
        self.interp_t2index = interp_t2index
        self.interp_p2dp = interp_p2dp
        self.interp_t2dt = interp_t2dt
        self.x2d = x2d
        self.y2d = y2d
        self.z2d = z2d


def get_1d_mesh_properties(x):
    """
    Return mesh properties based on an array of 1D mesh locations (assume monotonic).

    Assuming that the input x is the bounding exterior (main) mesh, this returns
    the interior (half) mesh locations and the cell centered spacings for each.

    The spacing BCs for cells at the boundary assume you are doing intergrals --> use half.

    :param x: 1D numpy array of grid locations.
    """
    # Compute the interior (half) mesh positions
    xh = 0.5*(x[1:] + x[0:-1])

    # Compute the spacing centered around the interior half mesh
    dxh = x[1:] - x[0:-1]

    # Compute the spacing centered on the bounding (main) mesh, interior first
    dx = xh[1:] - xh[0:-1]

    # Add the boundary points (half spacing)
    dx = np.concatenate([[xh[0] - x[0]], dx, [x[-1] - xh[-1]]])

    return xh, np.abs(dx), np.abs(dxh)

def downsamp_reg_grid(full_map, new_y, new_x, new_x_extents=None, new_y_extents=None,
                      full_x_extents=None, full_y_extents=None,
                      image_method=0, chd_method=0, periodic_x=True, y_units='sinlat', uniform_poles=True,
                      single_origin_image=None, uniform_no_data=True, out_full_latlon=True, sparse_wghts=False):
    """
    Input a PSI-map object (full_map) and re-sample to the new x, y coords.

    Grid definition: Assumes grids are regular, but non-uniform.  Input grid is defined by
    full_map.x and full_map.y.  Output grid is defined by new_x and new_y.  When extents are
    not explicitly specified, it is assumed that the min/max of each axis defines the
    extent of that axis (which implies a half-pixel at each edge).  All other axis values
    are assumed to be grid-centers.  When extents are specified, all axis values are assumed
    to be grid/pixel centers.

    :param full_map: PsiMap object
                The full-map (covers full sphere) to be downsampled
    :param new_y: Array like
                  Vector of pixel centers. Units can be specified with y_units.
                  Default is sin(lat).
    :param new_x: Array like
                  Vector of pixel centers. Currently only support longitude in
                  radians (phi).
    :param new_x_extents: list-like; len(new_x_extents) = 2
                          Depending on grid-definition, the x-extents of the output
                          grid may exceed the extents of 'new_x'. If None, assume
                          new_x_extents = [np.min(new_x), np.max(new_x)].
    :param new_y_extents: list-like; len(new_y_extents) = 2
                          Depending on grid-definition, the y-extents of the output
                          grid may exceed the extents of 'new_y'. If None, assume
                          new_x_extents = [np.min(new_y), np.max(new_y)].
    :param full_x_extents: list-like; len(full_x_extents) = 2
                           Depending on grid-definition, the x-extents of the input
                           grid may exceed the extents of full_map.x. If None, assume
                           full_x_extents = [np.min(full_map.x), np.max(full_map.x)].
    :param full_y_extents: list-like; len(full_y_extents) = 2
                           Depending on grid-definition, the y-extents of the input
                           grid may exceed the extents of full_map.y. If None, assume
                           full_y_extents = [np.min(full_map.y), np.max(full_map.y)].
    :param image_method: integer
                         0 - Average across overlapped pixels to downsample
                         1 - Use random sampling to downsample (not yet supported)
    :param chd_method: integer
                       0 - Average across overlapped pixels to downsample
    :param periodic_x: True/False
                       Should the left and right edge be treated as periodic
    :param y_units: character string
                    Latitude units: 'sinlat', 'lat_rad', 'lat_deg', 'theta'
    :param uniform_poles: True/False
                          Enforce uniformity at the north and south pole (across
                          the top and bottom edge of the map)
    :param single_origin_image: integer (default: None)
                                If the map contains information from only one image,
                                this value will be used to fill 'origin_image' in
                                the output map. If None, 'origin_image' will be set
                                to None.
    :param uniform_no_data: True/False
                            The algorithm is significantly sped-up if 'no_data' pixel
                            locations are consistent across all map grids (data, chd,
                            mu).  When 'no_data' pixels vary from grid to grid, set
                            to False.
    :param sparse_wghts: True/False
                         Use sparse matrices to hold/apply the reduction operators.
    :return: PsiMap object
             The new map object with grid defined by new_x and new_y.
    """
    # Determine/check the output grid extents
    if new_x_extents is not None:
        valid_extents_x = (np.min(new_x) >= new_x_extents[0]) & (np.max(new_x) <= new_x_extents[1])
        if not valid_extents_x:
            raise ValueError("The new grid defined by 'new_x' does not lie entirely within 'new_x_extents'")
        new_x_extents = np.array(new_x_extents)
    else:
        new_x_extents = np.array([np.min(new_x), np.max(new_x)])
    if new_y_extents is not None:
        valid_extents_y = (np.min(new_y) >= new_y_extents[0]) & (np.max(new_y) <= new_y_extents[1])
        if not valid_extents_y:
            raise ValueError("The new grid defined by 'new_y' does not lie entirely within 'new_y_extents'")
        new_y_extents = np.array(new_y_extents)
    else:
        new_y_extents = np.array([np.min(new_y), np.max(new_y)])

    # Determine/check the input grid extents
    if full_x_extents is not None:
        valid_extents_x = (np.min(full_map.x) >= full_x_extents[0]) & (np.max(full_map.x) <= full_x_extents[1])
        if not valid_extents_x:
            raise ValueError("The new grid defined by 'full_map.x' does not lie entirely within 'full_x_extents'")
        full_x_extents = np.array(full_x_extents)
    else:
        full_x_extents = np.array([np.min(full_map.x), np.max(full_map.x)])
    if full_y_extents is not None:
        valid_extents_y = (np.min(full_map.y) >= full_y_extents[0]) & (np.max(full_map.y) <= full_y_extents[1])
        if not valid_extents_y:
            raise ValueError("The new grid defined by 'full_map.y' does not lie entirely within 'full_y_extents'")
        full_y_extents = np.array(full_y_extents)
    else:
        full_y_extents = np.array([np.min(full_map.y), np.max(full_map.y)])

    # check that new coord range does not exceed old coord range
    x_in_range = (full_x_extents[0] <= new_x_extents[0]) & (full_x_extents[1] >= new_x_extents[1])
    y_in_range = (full_y_extents[0] <= new_y_extents[0]) & (full_y_extents[1] >= new_y_extents[1])
    if not x_in_range or not y_in_range:
        raise ValueError("The new grid extents exceed the"
                         "range of the existing grid extents.")

    start_time = time.time()
    # generate MapMesh object for both grids
    # new_theta = -np.arcsin(new_y) + np.pi/2
    # new_mesh = MapMesh(new_x, new_theta)
    if y_units == "sinlat":
        full_theta = -np.arcsin(full_map.y) + np.pi/2
        full_theta_extents = -np.arcsin(full_y_extents) + np.pi / 2
        sin_lat = full_map.y
        new_sin_lat = new_y
        new_sinlat_extents = new_y_extents
        new_theta = -np.arcsin(new_y) + np.pi/2
    elif y_units == "lat_rad":
        full_theta = -full_map.y + np.pi/2
        full_theta_extents = -full_y_extents + np.pi/2
        sin_lat = np.sin(full_map.y)
        new_sin_lat = np.sin(new_y)
        new_sinlat_extents = np.sin(new_y_extents)
        new_theta = -new_y + np.pi/2
    elif y_units == "lat_deg":
        full_theta = -(np.pi/180)*full_map.y + np.pi/2
        full_theta_extents = -(np.pi/180)*full_y_extents + np.pi/2
        sin_lat = np.sin((np.pi/180)*full_map.y)
        new_sin_lat = np.sin((np.pi/180)*new_y)
        new_sinlat_extents = np.sin((np.pi/180)*new_y_extents)
        new_theta = -(np.pi/180)*new_y + np.pi/2
    else:
        full_theta = full_map.y
        full_theta_extents = full_y_extents
        sin_lat = np.sin(-full_map.y + np.pi/2)
        new_sin_lat = np.sin(-new_y + np.pi/2)
        new_sinlat_extents = np.sin(-new_y_extents + np.pi/2)
        new_theta = new_y

    ## determine area of each CR map pixel
    phi_half = np.diff(full_map.x)/2 + full_map.x[:-1]
    lower_phi = np.append(full_x_extents[0], phi_half)
    upper_phi = np.append(phi_half, full_x_extents[1])
    diff_phi = upper_phi - lower_phi

    if y_units == "sinlat":
        sinlat_half = np.diff(full_map.y)/2 + full_map.y[:-1]
        theta_half = -np.arcsin(sinlat_half) + np.pi/2
    else:
        theta_half = np.diff(full_theta)/2 + full_theta[:-1]
    lower_theta = np.append(full_theta_extents[0], theta_half)
    upper_theta = np.append(theta_half, full_theta_extents[1])
    diff_theta = np.abs(upper_theta - lower_theta)
    full_sinlat_extents = np.sin(-full_theta_extents + np.pi/2)

    diff_phi_mat, diff_theta_mat = np.meshgrid(diff_phi, diff_theta)
    da = np.matmul(np.diag(np.sin(full_theta)), diff_theta_mat * diff_phi_mat)

    ## generate bin edges for each grid
    if y_units == "sinlat":
        new_y_interior_edges = (new_sin_lat[0:-1]+new_sin_lat[1:])/2
    else:
        new_y_interior_theta = (new_theta[0:-1] + new_theta[1:])/2
        new_y_interior_edges = np.sin(-new_y_interior_theta + np.pi/2)
    new_y_edges = np.concatenate([[new_sinlat_extents[0]], new_y_interior_edges, [new_sinlat_extents[1]]])
    new_x_interior_edges = (new_x[0:-1] + new_x[1:])/2
    new_x_edges = np.concatenate([[new_x_extents[0]], new_x_interior_edges, [new_x_extents[1]]])
    if y_units == "sinlat":
        old_y_interior_edges = sinlat_half
    else:
        old_y_interior_edges = np.sin(-theta_half + np.pi/2)
    old_y_edges = np.concatenate([[full_sinlat_extents[0]], old_y_interior_edges, [full_sinlat_extents[1]]])
    old_x_interior_edges = (full_map.x[0:-1] + full_map.x[1:])/2
    old_x_edges = np.concatenate([[full_map.x[0]], old_x_interior_edges, [full_map.x[-1]]])

    # determine overlap weights for each row and column of new grid
    #   include area-weighting in row associations
    new_y_n = len(new_sin_lat)
    old_y_n = len(sin_lat)
    old_y_widths = np.diff(old_y_edges)
    if sparse_wghts:
        row_weight_mat = sparse.lil_matrix((new_y_n, old_y_n), dtype=full_map.data.dtype)
        row_da_weight = sparse.lil_matrix((new_y_n, old_y_n), dtype=full_map.data.dtype)
    else:
        row_weight_mat = np.zeros((new_y_n, old_y_n), dtype=full_map.data.dtype)
        row_da_weight = np.zeros((new_y_n, old_y_n), dtype=full_map.data.dtype)

    for new_y_index in range(new_y_n):
        # determine linear row-weighting of original pixels to new pixels
        # new_hist = np.ones(1)
        temp_edges = new_y_edges[new_y_index:(new_y_index+2)]
        # old_hist = hist_integration(new_hist, temp_edges, old_y_edges)
        pixel_portions = pixel_portion_overlap1D(temp_edges, old_y_edges)
        bin_indices = np.where(pixel_portions > 0.)
        bin_weights = pixel_portions[bin_indices]
        row_da_weight[new_y_index, bin_indices] = bin_weights
        # also weight by pixel width(height).
        # area_weights = da[bin_indices, 1]
        area_weights = old_y_widths[bin_indices]
        bin_weights = bin_weights*area_weights
        # normalize
        bin_weights = bin_weights/bin_weights.sum()
        # store indices and weights for each row
        row_weight_mat[new_y_index, bin_indices] = bin_weights

    # repeat for columns
    new_x_n = len(new_x)
    old_x_n = len(full_map.x)
    old_x_widths = np.diff(old_x_edges)
    if sparse_wghts:
        column_weight_mat = sparse.lil_matrix((old_x_n, new_x_n), dtype=full_map.data.dtype)
        col_da_weight = sparse.lil_matrix((old_x_n, new_x_n), dtype=full_map.data.dtype)
    else:
        column_weight_mat = np.zeros((old_x_n, new_x_n), dtype=full_map.data.dtype)
        col_da_weight = np.zeros((old_x_n, new_x_n), dtype=full_map.data.dtype)
    for new_x_index in range(new_x_n):
        # determine linear row-weighting of original pixels to new pixels
        # new_hist = np.ones(1)
        temp_edges = new_x_edges[new_x_index:(new_x_index + 2)]
        # old_hist = hist_integration(new_hist, temp_edges, old_x_edges)
        pixel_portions = pixel_portion_overlap1D(temp_edges, old_x_edges)
        bin_indices = np.where(pixel_portions > 0.)
        bin_weights = pixel_portions[bin_indices]
        col_da_weight[bin_indices, new_x_index] = bin_weights
        # multiply by pixel widths
        bin_weights = bin_weights * old_x_widths[bin_indices]
        # normalize
        bin_weights = bin_weights/bin_weights.sum()
        # store indices and weights for each column
        column_weight_mat[bin_indices, new_x_index] = bin_weights

    # prepare (de-linked) data for weighted-averaging
    full_data = full_map.data.copy()
    no_data_index = full_data == full_map.no_data_val
    full_data[no_data_index] = 0.
    # apply the row and column reduction by matrix multiplication
    if sparse_wghts:
        row_weight_mat = row_weight_mat.tocsr()
        column_weight_mat = column_weight_mat.tocsc()
        reduced_data = row_weight_mat @ full_data @ column_weight_mat
        # also calculate da in the new grid
        reduced_grid_da = row_da_weight @ da @ col_da_weight
    else:
        row_reduced_data = np.matmul(row_weight_mat, full_data)
        reduced_data = np.matmul(row_reduced_data, column_weight_mat)
        # also calculate da in the new grid
        reduced_grid_da = np.matmul(np.matmul(row_da_weight, da), col_da_weight)

    no_data_da = da.copy()
    no_data_da[no_data_index] = 0.
    if sparse_wghts:
        row_da_weight = row_da_weight.tocsr()
        col_da_weight = col_da_weight.tocsc()
        reduced_no_data_da = row_da_weight @ no_data_da @ col_da_weight
    else:
        reduced_no_data_da = np.matmul(np.matmul(row_da_weight, no_data_da),
                                       col_da_weight)
    # use the area ratio to improve intensity estimate at data boundaries (and
    # better estimate the boundary)
    da_ratio = reduced_no_data_da/reduced_grid_da
    new_no_data_index = da_ratio < 0.5
    reduced_data[new_no_data_index] = full_map.no_data_val
    reduced_data[~new_no_data_index] = reduced_data[~new_no_data_index]/ \
        da_ratio[~new_no_data_index]

    if single_origin_image is not None:
        origin_image = np.zeros(reduced_data.shape)
        origin_image[~new_no_data_index] = single_origin_image
    else:
        origin_image = None

    # now apply reduction to chd map
    if hasattr(full_map, "chd") and full_map.chd is not None:
        chd_data = full_map.chd.copy()
        if not uniform_no_data:
            # recalculate the no_data locations and area ratio
            no_chd_index = chd_data == full_map.no_data_val
            reduced_grid_da = np.matmul(np.matmul(row_da_weight, da), col_da_weight)
            no_data_da = da.copy()
            no_data_da[no_chd_index] = 0.
            reduced_no_data_da = np.matmul(np.matmul(row_da_weight, no_data_da),
                                           col_da_weight)
            # use the area ratio to improve intensity estimate at data boundaries (and
            # better estimate the boundary)
            da_ratio = reduced_no_data_da/reduced_grid_da
            new_no_chd_index = da_ratio < 0.5
        else:
            no_chd_index = no_data_index
            new_no_chd_index = new_no_data_index

        chd_data[no_chd_index] = 0.
        # apply the row and column reduction by matrix multiplication
        row_reduced_chd = np.matmul(row_weight_mat, chd_data)
        reduced_chd = np.matmul(row_reduced_chd, column_weight_mat)
        # use the area ratio to improve chd estimate at data boundaries (and
        # better estimate the boundary)
        reduced_chd[new_no_chd_index] = full_map.no_data_val
        reduced_chd[~new_no_chd_index] = reduced_chd[~new_no_chd_index] / \
            da_ratio[~new_no_chd_index]
    else:
        reduced_chd = None

    # now apply reduction to mu values
    if full_map.mu is not None:
        mu_data = full_map.mu.copy()
        # if not uniform_no_data:
        #     # recalculate the no_data locations and area ratio
        #     no_data_index = chd_data == full_map.no_data_val
        #     reduced_grid_da = np.matmul(np.matmul(row_da_weight, da), col_da_weight)
        #     no_data_da = da.copy()
        #     no_data_da[no_data_index] = 0.
        #     reduced_no_data_da = np.matmul(np.matmul(row_da_weight, no_data_da),
        #                                    col_da_weight)
        #     # use the area ratio to improve intensity estimate at data boundaries (and
        #     # better estimate the boundary)
        #     da_ratio = reduced_no_data_da/reduced_grid_da
        #     new_no_data_index = da_ratio < 0.5
        #
        # mu_data[no_data_index] = 0.
        # apply the row and column reduction by matrix multiplication
        row_reduced_mu = np.matmul(row_weight_mat, mu_data)
        reduced_mu = np.matmul(row_reduced_mu, column_weight_mat)
        # use the area ratio to improve mu estimate at data boundaries (and
        # better estimate the boundary)
        # reduced_mu[new_no_data_index] = full_map.no_data_val
        # reduced_mu[~new_no_data_index] = reduced_mu[~new_no_data_index]/ \
        #     da_ratio[~new_no_data_index]
    else:
        reduced_mu = None

    # now apply reduction to map_lon values (?)
    if hasattr(full_map, "map_lon") and full_map.map_lon is not None:
        map_lon = None
    else:
        map_lon = None

    if uniform_poles:
        no_data_vec = reduced_data[0, ] == full_map.no_data_val
        if np.any(~no_data_vec):
            reduced_data[0, ] = np.mean(reduced_data[0, ~no_data_vec])
        no_data_vec = reduced_data[-1, ] == full_map.no_data_val
        if np.any(~no_data_vec):
            reduced_data[-1, ] = np.mean(reduced_data[-1, ~no_data_vec])
        if reduced_chd is not None:
            no_data_vec = reduced_chd[0, ] == full_map.no_data_val
            if np.any(~no_data_vec):
                reduced_chd[0, ] = np.mean(reduced_chd[0, ~no_data_vec])
            no_data_vec = reduced_chd[-1, ] == full_map.no_data_val
            if np.any(~no_data_vec):
                reduced_chd[-1, ] = np.mean(reduced_chd[-1, ~no_data_vec])
        if reduced_mu is not None:
            no_data_vec = reduced_mu[0, ] == full_map.no_data_val
            if np.any(~no_data_vec):
                reduced_mu[0, ] = np.mean(reduced_mu[0, ~no_data_vec])
            no_data_vec = reduced_mu[-1, ] == full_map.no_data_val
            if np.any(~no_data_vec):
                reduced_mu[-1, ] = np.mean(reduced_mu[-1, ~no_data_vec])

    if periodic_x:
        new_x_widths = np.diff(new_x_edges)
        # average half-pixels from left and right edge
        reduced_data[:, 0], reduced_data[:, -1] = periodic_x_avg(
            reduced_data[:, 0], reduced_data[:, -1], new_x_widths, full_map.no_data_val)

        if reduced_mu is not None:
            reduced_mu[:, 0], reduced_mu[:, -1] = periodic_x_avg(
                reduced_mu[:, 0], reduced_mu[:, -1], new_x_widths, full_map.no_data_val)

        if reduced_chd is not None:
            reduced_chd[:, 0], reduced_chd[:, -1] = periodic_x_avg(
                reduced_chd[:, 0], reduced_chd[:, -1], new_x_widths, full_map.no_data_val)

    end_time = time.time()
    # print(end_time - start_time, " seconds elapsed.\n")

    # quick plot testing (remove at clean-up)
    # import matplotlib.pyplot as plt
    # plt.figure(0)
    # plt.imshow(full_map.data, origin='lower')
    # plt.title("Original Map")
    # plt.figure(1)
    # plt.imshow(reduced_data, origin='lower')
    # plt.title("Reduced Map")

    if full_map.method_info is not None:
        # alter method 'GridSize_sinLat' to new resolution
        method_info = full_map.method_info.copy()
        y_index = method_info.meth_name.eq('GridSize_sinLat') & \
            method_info.var_name.eq('n_SinLat')
        method_info.loc[y_index, 'var_val'] = new_y_n
        x_index = method_info.meth_name.eq('GridSize_sinLat') & \
            method_info.var_name.eq('n_phi')
        method_info.loc[x_index, 'var_val'] = new_x_n
    else:
        method_info = None

    # generate new map object and fill
    map_type = type(full_map).__name__
    if map_type == "PsiMap":
        new_map = psi_d_types.PsiMap(data=reduced_data, x=new_x, y=new_y, mu=reduced_mu,
                                     origin_image=origin_image, map_lon=map_lon, chd=reduced_chd,
                                     no_data_val=full_map.no_data_val)
        new_map.append_method_info(method_info)
        new_map.append_map_info(full_map.map_info)
        new_map.append_data_info(full_map.data_info)
    elif map_type == "MagnetoMap":
        new_map = psi_d_types.MagnetoMap(data=reduced_data, x=new_x, y=new_y, mu=reduced_mu,
                                         no_data_val=full_map.no_data_val)
        # may need to add additional database info pass-through here
        new_map.sunpy_meta = full_map.sunpy_meta
    else:
        new_map = None

    return new_map


def pixel_portion_overlap1D(edges, new_edges):
    # function to calc new pixel overlap with a single pixel
    """

    :param edges: list-like with two entries (sorted increasing)
           Original pixel edges
    :param new_edges: list-like (sorted increasing)
           New pixel edges
    :return: np.ndarray vector with length=len(new_edges)-1
             The portion of each new pixel that overlaps the original pixel.
    """
    # initiate results vector
    n_bins = len(new_edges) - 1
    out_vec = np.zeros(n_bins)

    left_edges = new_edges[0:-1]
    right_edges = new_edges[1:]
    left_in = left_edges < edges[1]
    left_out = left_edges < edges[0]
    right_in = right_edges > edges[0]
    right_out = right_edges > edges[1]

    # for extremely large arrays (approx n_bins>2E5), this is probably faster
    # temp_index = np.searchsorted(left_edges, edges[1], side='right')
    # left_in = np.zeros(n_bins, dtype=bool)
    # left_in[:temp_index] = True

    consume_bin = left_out & right_out
    # check for single new pixel that surrounds old pixel
    if any(consume_bin):
        # calculate portion of surrounding bin that overlaps old pixel
        out_vec[consume_bin] = (edges[1] - edges[0])/(right_edges[consume_bin] -
                                                      left_edges[consume_bin])
        return out_vec

    # check left overlap for partial overlap
    left_overlap = left_out & right_in
    out_vec[left_overlap] = (right_edges[left_overlap] - edges[0]) / \
                            (right_edges[left_overlap] - left_edges[left_overlap])

    # check for partial right overlap
    right_overlap = right_out & left_in
    out_vec[right_overlap] = (edges[1] - left_edges[right_overlap]) / \
                             (right_edges[right_overlap] - left_edges[right_overlap])

    # remaining overlap pixels fall inside original pixel
    full_overlap = ~left_out & ~right_out
    out_vec[full_overlap] = 1.

    return out_vec

def periodic_x_avg(l_vec, r_vec, x_edges, no_data_val):

    # average half-pixels from left and right edge
    left_no_data = l_vec == no_data_val
    right_no_data = r_vec == no_data_val

    both_index = ~left_no_data & ~right_no_data
    average_vec = (x_edges[0]*l_vec[both_index] +
                   x_edges[-1]*r_vec[both_index])/ \
                  (x_edges[0] + x_edges[-1])
    l_vec[both_index] = average_vec
    r_vec[both_index] = average_vec

    right_index = left_no_data & ~right_no_data
    l_vec[right_index] = r_vec[right_index]
    left_index = ~left_no_data & right_no_data
    r_vec[left_index] = l_vec[left_index]

    return l_vec, r_vec
