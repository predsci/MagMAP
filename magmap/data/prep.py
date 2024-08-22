"""
Modules for manipulating raw images (fits) to our to prep for mapping.
- Here the relevant data reduction and preparation steps are applied, including:
  - Image rotation
  - Use FITS header to generate axis scales
"""
import numpy as np
import astropy.units as u


def rotate_map_nopad(raw_map):
    """
    Wrapper for sunpy.map.mapbase.rotate that does rotation how we like it.
    - The images are recentered and the extra padding produced by rotate is
      removed after rotation (got this from the aia_prep routine in Sunpy 1.03).

    :param raw_map: Raw map from FITS file
    :return newmap: New map rotated to polar-north=up and padding removed
    """
    tempmap = raw_map.rotate(recenter=True)

    # extract center from padded map.rotate output
    # - crpix1 and crpix2 will be equal (recenter=True) -> does not work with submaps
    center = tempmap.meta['crpix1']
    # newer sunpy wants bottom_left and top_right rather than axis1_range and axis2_range
    if (center % 1.) == 0:
        # Implies an odd number of pixels, remove an extra pixel from top and right
        # Assumes original shape is even
        bottom_left = (center - np.array([1, 1])*raw_map.data.shape[0]/2 +
                       np.array([0, 0]))*u.pix
        top_right = (center + np.array([1, 1])*raw_map.data.shape[0]/2 -
                     np.array([1, 1]))*u.pix
    else:
        bottom_left = (center - np.array([1, 1])*raw_map.data.shape[0]/2)*u.pix
        top_right = (center + np.array([1, 1])*raw_map.data.shape[0]/2)*u.pix

    newmap = tempmap.submap(bottom_left=bottom_left, top_right=top_right)

    return newmap


def get_scales_from_fits(fits_meta):
    """
    Compute the solar X and solar Y 1D scale arrays from the fits metadata.

    Return scales in solar radii and assume that image-up and image-right are positive.
    :param fits_meta:   sunpy.util.metadata.MetaDict
                        Meta data loaded from fits file by example = sunpy.map.Map()
                        Meta data is accessed by 'example.meta'
    :return:    tuple of arrays
                First array is x-axis of image space in solar radii.
                Second array is y-axis of image space in solar radii.
    """

    if 'rsun_arc' in fits_meta.keys():
        arcsec_per_radii = fits_meta['rsun_arc']
    else:
        arcsec_per_radii = fits_meta['rsun_obs']
    # y-axis pars
    crpix2 = fits_meta['crpix2']
    cdelt2 = fits_meta['cdelt2']
    naxis2 = fits_meta['naxis2']
    # x-axis pars
    crpix1 = fits_meta['crpix1']
    cdelt1 = fits_meta['cdelt1']
    naxis1 = fits_meta['naxis1']

    # pixel locations (starting at 1 and not 0, per the fits standard)
    xpix_num = np.arange(start=1, stop=naxis1+1, step=1)
    rel_xpix = xpix_num - crpix1
    # convert to arcsec
    x_arcsec = rel_xpix*cdelt1
    # convert to solar radii
    x_radii = x_arcsec/arcsec_per_radii

    # pixel locations (starting at 1 and not 0, per the fits standard)
    ypix_num = np.arange(start=1, stop=naxis2+1, step=1)
    rel_ypix = ypix_num - crpix2
    # convert to arcsec
    y_arcsec = rel_ypix * cdelt2
    # convert to solar radii
    y_radii = y_arcsec / arcsec_per_radii

    return x_radii, y_radii


def get_scales_from_fits_map(fits_meta):
    """
    Compute the solar X and solar Y 1D scale arrays from the fits metadata.

    Written specifically for HMI_Mrmap_latlon_720s fits files, but should generally
    work to return X/Y coordinates in the native units.
    :param fits_meta:   sunpy.util.metadata.MetaDict
                        Meta data loaded from fits file by example = sunpy.map.Map()
                        Meta data is accessed by 'example.meta'
    :return:    tuple of arrays
                First array is x-axis of image space in native units.
                Second array is y-axis of image space in native units.
    """

    # y-axis pars
    crpix2 = fits_meta['crpix2']
    cdelt2 = fits_meta['cdelt2']
    naxis2 = fits_meta['naxis2']
    crval2 = fits_meta['crval2']
    # x-axis pars
    crpix1 = fits_meta['crpix1']
    cdelt1 = fits_meta['cdelt1']
    naxis1 = fits_meta['naxis1']
    crval1 = fits_meta['crval1']

    # pixel locations (starting at 1 and not 0, per the fits standard)
    xpix_num = np.arange(start=1, stop=naxis1+1, step=1)
    rel_xpix = xpix_num - crpix1
    # convert to native units
    x_native = rel_xpix*cdelt1 + crval1

    # pixel locations (starting at 1 and not 0, per the fits standard)
    ypix_num = np.arange(start=1, stop=naxis2+1, step=1)
    rel_ypix = ypix_num - crpix2
    # convert to native units
    y_native = rel_ypix * cdelt2 + crval2

    return x_native, y_native

