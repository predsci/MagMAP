"""
Set of functions to change coordinates and interpolate between images and maps
"""

import numpy as np
import scipy.interpolate as sp_interp

import oftpy.utilities.datatypes.datatypes as psi_dt
# import astropy_healpix


def c2s(x, y, z):
    """
    convert numpy arrays of x,y,z (cartesian) to r,t,p (spherical)
    """
    x2 = x**2
    y2 = y**2
    z2 = z**2
    r = np.sqrt(x2 + y2 + z2)
    t = np.arctan2(np.sqrt(x2 + y2), z)
    p = np.arctan2(y, x)

    # arctan2 returns values from -pi to pi but I want 0-2pi --> use fmod
    twopi = np.pi*2
    p = np.fmod(p + twopi, twopi)

    return r, t, p


def s2c(r, t, p):
    """
    convert numpy arrays of r,t,p (spherical) to x,y,z (cartesian)
    """
    ct = np.cos(t)
    st = np.sin(t)
    cp = np.cos(p)
    sp = np.sin(p)
    x = r*cp*st
    y = r*sp*st
    z = r*ct
    return x, y, z


def get_arclength(chord_length, radius=1.0):
    """
    convert the length of the chord connecting two points on a great circle to the
    arc length connecting them on the great circle.
    """
    arclength = 2*radius*np.arcsin(0.5*chord_length/radius)
    return arclength


def map_grid_to_image(map_x, map_y, R0=1.0, obsv_lon=0.0, obsv_lat=0.0, image_crota2=0.0):
    """
    Given a set of xy coordinate pairs, in map-space (x:horizontal phi axis, y:vertical sin(theta) axis, radius:R0),
    rotate and change variables to image space (x:horizontal, y:vertical, z:coming out of image, theta=0 at center of
    image, phi=0 at x=R0)
    :param map_x: numpy 1D array of map pixel x-locations [0, 2pi]
    :param map_y: numpy 1D array of map pixel y-locations [-1, 1]
    :param R0: Image radius in solar radii
    :param obsv_lon: Carrington longitude of image observer
    :param obsv_lat: Carrington latitude of image observer
    :param image_crota2: Degrees counterclockwise rotation needed to put solar-north-up in the image. This is generally
    a parameter in the .fits metadata named 'crota2'.
    :return:
    """

    # get map 3D spherical coords
    map_theta = np.pi / 2 - np.arcsin(map_y)
    map_phi = map_x

    # get map 3D cartesian coords
    map3D_x = R0 * np.sin(map_theta) * np.cos(map_phi)
    map3D_y = R0 * np.sin(map_theta) * np.sin(map_phi)
    map3D_z = R0 * np.cos(map_theta)

    # generate rotation matrix (from Carrington to image solar-north-up)
    rot_mat1 = map_to_image_rot_mat(obsv_lon, obsv_lat)
    # generate rotation matrix (from image solar-north-up to image orientation)
    rot_mat2 = snu_to_image_rot_mat(image_crota2)
    # combine rotations
    rot_mat = np.matmul(rot_mat2, rot_mat1)
    # construct coordinate array
    coord_array = np.array([map3D_x, map3D_y, map3D_z])
    # apply rotation matrix to coordinates
    image3D_coord = np.matmul(rot_mat, coord_array)

    # numeric error occasionally results in |z| > R0. This is a problem for np.arccos()
    image3D_coord[2, image3D_coord[2, :] > R0] = R0
    image3D_coord[2, image3D_coord[2, :] < -R0] = -R0
    # a more proper solution is to re-normalize each coordinate set, but that would be
    # more expensive (to fix error on an order we are not worried about). Also, the
    # numeric error of the renormalization could still create problems with arccos().

    image_phi = np.arctan2(image3D_coord[1, :], image3D_coord[0, :])
    image_theta = np.arccos(image3D_coord[2, :] / R0)

    return image3D_coord[0, :], image3D_coord[1, :], image3D_coord[2, :], image_theta, image_phi


def image_grid_to_CR(image_x, image_y, R0=1.0, obsv_lat=0, obsv_lon=0, get_mu=False, outside_map_val=-9999.):
    """
    Given vector coordinate pairs in solar radii units and the observer angles, transform to map coords.
    :param image_x: vector of x coordinates
    :param image_y: vector of y coordinates
    :param R0: Assumed radius in solar radii.
    :param obsv_lat: Carrington latitude (degrees from equator) of observing instrument
    :param obsv_lon: Carrington longitude (degrees) of observing instrument
    :return:
    """

    # for images, we assume that the z-axis is perpendicular to the image plane, in the direction
    # of the observer, and located at the center of the image.

    # mask points outside of R0
    use_index = image_x ** 2 + image_y ** 2 <= R0 ** 2
    use_x = image_x[use_index]
    use_y = image_y[use_index]

    # Find z coord (we can assume it is in the positive direction)
    # use_z = np.sqrt(R0 ** 2 - use_x ** 2 - use_y ** 2)
    # to be numerically equivalent to the use_index definition, change to this:
    use_z = np.sqrt(R0 ** 2 - (use_x ** 2 + use_y ** 2))

    # Calc image_theta, image_phi, and image_mu
    if get_mu:
        image_mu = np.full(image_x.shape, outside_map_val)
        # image_phi = np.arctan2(image_y, image_x)
        use_theta = np.arccos(use_z / R0)
        image_mu[use_index] = np.cos(use_theta)

    # generate map-to-image rotation matrix
    rot_mat = map_to_image_rot_mat(obsv_lon, obsv_lat)
    # invert/transpose for image-to-map rotation matrix
    rev_rot = rot_mat.transpose()
    # construct coordinate array
    coord_array = np.array([use_x, use_y, use_z])
    # apply rotation matrix to coordinates
    map3D_coord = np.matmul(rev_rot, coord_array)

    # Occasionally numeric error from the rotation causes a z magnitude to be greater than R0
    num_err_z_index = np.abs(map3D_coord[2, :]) > R0
    map3D_coord[2, num_err_z_index] = np.sign(map3D_coord[2, num_err_z_index]) * R0
    # Convert map cartesian to map theta and phi
    cr_theta = np.arccos(map3D_coord[2, :] / R0)
    cr_phi = np.arctan2(map3D_coord[1, :], map3D_coord[0, :])
    # Change phi range from [-pi,pi] to [0,2pi]
    neg_phi = cr_phi < 0
    cr_phi[neg_phi] = cr_phi[neg_phi] + 2 * np.pi

    cr_theta_all = np.full(image_x.shape, outside_map_val)
    cr_phi_all = np.full(image_x.shape, outside_map_val)

    cr_theta_all[use_index] = cr_theta
    cr_phi_all[use_index] = cr_phi

    if get_mu:
        return cr_theta_all, cr_phi_all, image_mu
    else:
        return cr_theta_all, cr_phi_all, None


def interpolate2D_regular2irregular(reg_x, reg_y, reg_vals, eval_x, eval_y):
    """
    For a 2D MxN regular grid, interpolate values to the K grid points in eval_x and eval_y.
    :param reg_x: numpy vector of x-coordinates (length N)
    :param reg_y: numpy vector of y-coordinates (length M)
    :param reg_vals: numpy MxN array_like containing the grid values
    :param eval_x: numpy column vector (length K) of x-coordinates to evaluate at.
    :param eval_y: numpy column vector (length K) of y-coordinates to evaluate at.
    :return: vector length K of interpolation results.
    """
    # Setup interpolation function and grd to evaluate on
    interp_fn = sp_interp.RegularGridInterpolator((reg_x, reg_y), reg_vals)
    eval_pts = np.array([eval_y, eval_x]).transpose()

    # Do interpolation
    interp_result_vec = interp_fn(eval_pts)

    return interp_result_vec


def interp_los_image_to_map(image_in, R0, map_x, map_y, no_data_val=-9999.):
    map_nxcoord = len(map_x)
    map_nycoord = len(map_y)

    # initialize grid to receive interpolation with values of NaN
    interp_result = np.full((map_nycoord, map_nxcoord), no_data_val, dtype='<f4')

    # convert 1D map axis to full list of coordinates
    mat_x, mat_y = np.meshgrid(map_x, map_y)
    # convert matrix of coords to vector of coords (explicitly select row-major vectorizing)
    map_x_vec = mat_x.flatten(order="C")
    map_y_vec = mat_y.flatten(order="C")
    interp_result_vec = interp_result.flatten(order="C")
    # determine if image is solar-north-up, or needs an additional rotation
    if hasattr(image_in, "sunpy_meta"):
        if "crota2" in image_in.sunpy_meta.keys():
            image_crota2 = image_in.sunpy_meta['crota2']
        else:
            image_crota2 = 0.
    else:
        image_crota2 = 0.
    # convert map grid variables to image space
    image_x, image_y, image_z, image_theta, image_phi = map_grid_to_image(map_x_vec, map_y_vec, R0=R0,
                                                                          obsv_lon=image_in.info['cr_lon'],
                                                                          obsv_lat=image_in.info['cr_lat'],
                                                                          image_crota2=image_crota2)
    # only interpolate points on the front half of the sphere
    interp_index = image_z > 0

    if type(image_in) is psi_dt.IITImage:
        im_data = image_in.iit_data
    elif type(image_in) is psi_dt.LBCCImage:
        im_data = image_in.lbcc_data
    else:
        im_data = image_in.data
    interp_vec = interpolate2D_regular2irregular(image_in.x, image_in.y, im_data, image_x[interp_index],
                                                 image_y[interp_index])
    interp_result_vec[interp_index] = interp_vec
    # reformat result to matrix form
    interp_result = interp_result_vec.reshape((map_nycoord, map_nxcoord), order="C")

    mu_vec = np.cos(image_theta)
    mu_mat = mu_vec.reshape((map_nycoord, map_nxcoord), order="C")

    out_obj = psi_dt.InterpResult(interp_result, map_x, map_y, mu_mat=mu_mat)

    return out_obj


def map_to_image_rot_mat(obsv_lon, obsv_lat):
    # del_phi = -obsv_lon*np.pi/180 - np.pi/2
    # int3D_x = np.cos(del_phi)*map3D_x - np.sin(del_phi)*map3D_y
    # int3D_y = np.cos(del_phi)*map3D_y + np.sin(del_phi)*map3D_x
    # int3D_z = map3D_z

    # rotate phi (about map z-axis. observer phi goes to -y)
    del_phi = -obsv_lon * np.pi / 180 - np.pi / 2
    rot1 = np.array([[np.cos(del_phi), -np.sin(del_phi), 0.],
                     [np.sin(del_phi), np.cos(del_phi), 0.], [0., 0., 1.], ])

    # rotate theta (about x-axis. observer theta goes to +z)
    del_theta = obsv_lat * np.pi / 180 - np.pi / 2
    rot2 = np.array([[1., 0., 0.], [0., np.cos(del_theta), -np.sin(del_theta)],
                     [0., np.sin(del_theta), np.cos(del_theta)]])

    tot_rot = np.matmul(rot2, rot1)

    return tot_rot

def snu_to_image_rot_mat(crota2):
    # Use the 'crota2' parameter (degrees counterclockwise) from fits metadata to rotate points
    # from solar-north-up (snu) to image orientation.
    # Assumes that we are in 3-D image-coordinates and rotating about
    # the observer line-of-sight (image z-axis)
    # Also assumes that the image-space horizontal axis is 'x' and increases to the right and
    # that the image-space vertical axis is 'y' and increases up.  Positive z-axis is toward
    # the observer.

    # convert to radians
    crota_rad = np.pi*crota2/180
    rot_mat = np.array([[np.cos(crota_rad), -np.sin(crota_rad), 0.],
                        [np.sin(crota_rad), np.cos(crota_rad), 0.],
                        [0., 0., 1.]])
    return rot_mat

