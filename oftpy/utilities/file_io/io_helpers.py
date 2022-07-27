"""
Generic helper routines.

Includes routines for downloading, organizing, and manipulating FITS files.

These might be split up into different modules later.
"""
import re

from six.moves.urllib.request import urlretrieve
from six.moves.urllib.error import HTTPError, URLError
import os
import pandas as pd
from collections import OrderedDict
import astropy.io.fits
import astropy.time as astro_time
import numpy as np
import sunpy.coordinates.sun


def download_url(url, fpath, overwrite=False, verbose=True):
    """
    Quick function to download a file specified by a url. I used drms/client.py as am
    example to strip out the exact functionality that I needed.
    - This is NOT very general --> don't parse the url for a path/file, just require that
      the call specifies one.
    - If file downloaded, this returns 0
    - If download error, this returns 1
    - If file already exists (no overwrite), this returns 2
    - If file already exists (overwrite), this returns 3
    """

    exit_flag = 0
    # check if the file exists first
    if os.path.isfile(fpath):
        if not overwrite:
            print("    " + fpath + " exists! SKIPPING!")
            return 2
        else:
            print("    " + fpath + " exists! OVERWRITING!")
            exit_flag = 3

    # use a temporary filename during download
    fpath_tmp = fpath + '.part'
    try:
        if verbose:
            print('### Downloading file:')
            print('  url:  ' + url)
            print('  path: ' + fpath)
        urlretrieve(url, fpath_tmp)
    except (HTTPError, URLError):
        print('  -> Error: Could not download file')
        return 1
    except ConnectionResetError:
        print('  -> Error: Connection was reset')
        return 1
    else:
        os.rename(fpath_tmp, fpath)
        if verbose:
            print('  Done!')
        return exit_flag


def construct_path_and_fname(base_dir, dtime, prefix, postfix, extension, inst=None, mkdir=True):
    """
    Quick function to build a subdirectory path and filename for saving/reading data
    - This parses a builtin datetime object to get the time strings used to build the info
    - the prefix and postfix should typically be instrument/series name and the wavelength respectively
    - it returns the subdirectory path and filename
    """

    # parse the datetime object:
    YYYY = '{:0>4}'.format(str(dtime.year))
    MM = '{:0>2}'.format(str(dtime.month))
    DD = '{:0>2}'.format(str(dtime.day))
    HH = '{:0>2}'.format(str(dtime.hour))
    NN = '{:0>2}'.format(str(dtime.minute))
    SS = '{:0>2}'.format(str(dtime.second))

    # build the subdirectory path
    sub_dir = os.path.join(base_dir, YYYY, MM, DD)

    # make the directory if needed
    if mkdir:
        # first check if the main directory exists
        if not os.path.isdir(base_dir):
            raise Exception('Base path does not exist! ' + base_dir)
            return None, None
        # check if the subdirectory exists
        if not os.path.isdir(sub_dir):
            os.makedirs(sub_dir)

    # build the filename
    if inst is not None:
        fname = prefix + '_' + inst + "_" + YYYY + MM + DD + 'T' + HH + NN + SS + '_' + postfix + '.' + extension
    elif postfix == "":
        fname = prefix + '_' + YYYY + MM + DD + 'T' + HH + NN + SS + '.' + extension
    else:
        fname = prefix + '_' + YYYY + MM + DD + 'T' + HH + NN + SS + '_' + postfix + '.' + extension

    return sub_dir, fname


def construct_hdf5_pre_and_post(chd_meta):
    """
    Standardize/Centralize hdf5 image filename production
    :param chd_meta: image meta dictionary. Output of euv_utils.py - get_metadata()
    :return: prefix, postfix, and extension strings
    """
    prefix = chd_meta['instrument'].lower().replace('-', '') + '_lvl2'
    postfix = str(chd_meta['wavelength'])
    extension = 'h5'

    return prefix, postfix, extension


def custom_dataframe(times, jds, urls, spacecraft, instrument, filters):
    """
    General function designed to take information from any query and turn it into a sliceable
    pandas dataframe with only the information I want for sorting/downloading
    The basic idea here is to make it easier to work with AIA and EUVI query results
    :return:
    """
    data_frame = pd.DataFrame(OrderedDict({'spacecraft': spacecraft, 'instrument': instrument,
                                           'filter': filters, 'time': times, 'jd': jds, 'url': urls}))

    return data_frame


def compress_uncompressed_fits_image(infile, outfile):
    """
    read an uncompressed fits file and compress it using the default rice compression
    - you can either save it to a different file or overwrite it by supplying the same fname
    - The silentfix is to avoid warnings/and or crashes when astropy encounters nans in the header values
    """
    hdulist = astropy.io.fits.open(infile)
    hdulist.verify('silentfix')
    hdr = hdulist[0].header

    # write out the file
    comp_hdu = astropy.io.fits.CompImageHDU(hdulist[0].data, hdr)
    hdulist.close()
    comp_hdu.writeto(outfile, output_verify='silentfix', overwrite=True, checksum=True)


def uncompress_compressed_fits_image(infile, outfile, int=False):
    """
    read an compressed fits file and uncompress it so that any .fits reader can read it.
    - you can either save it to a different file or overwrite it by supplying the same fname
    - The silentfix is to avoid warnings/and or crashes when astropy encounters nans in the header values
    """
    hdulist = astropy.io.fits.open(infile)
    hdulist.verify('silentfix')

    # check that it is a compressed image by looking at the length of the hdulist
    if len(hdulist) != 2:
        hdulist.info()
        raise Exception("This FITS file does not look like a simple compressed image!")

    hdr = hdulist[1].header

    # By default Astropy will convert the compressed data to float64
    data = hdulist[1].data

    # for some files (e.g. un-prepped STEREO you might want unsigned integer)
    if int:
        data = data.astype(np.uint16)

    # write out the file
    hdu = astropy.io.fits.PrimaryHDU(data, hdr)
    hdulist.close()
    hdu.writeto(outfile, output_verify='silentfix', overwrite=True, checksum=True)


def read_uncompressed_fits_image(infile):
    """
    read an uncompressed fits file and return the data and header
    - The silentfix is to avoid warnings/and or crashes when astropy encounters nans in the header values
    """
    hdulist = astropy.io.fits.open(infile)
    hdulist.verify('silentfix')
    hdr = hdulist[0].header
    data = hdulist[0].data

    return data, hdr


def write_sunpy_map_as_fits(outfile, map, dtype=np.uint16):
    """
    Helper function to take a sunpy map object and save the data as a fits file.

    - This function circumvents sunpy.io in order to gain flexibility in the output
      data type and how the fits library is called.

    - We use it mainly to create a file that looks like raw STEREO EUVI
      images in unsigned int format
    """
    # get the numpy array of image data, convert it to the desired dtype.
    data = dtype(map.data)

    # get the fits header (an astropy Header object)
    header = map.fits_header

    # start a new fits file object
    hdu = astropy.io.fits.PrimaryHDU(data, header=header)

    # build the hdu list
    hdulist = astropy.io.fits.HDUList([hdu])

    # write the file
    hdulist.close()
    hdu.writeto(outfile, output_verify='silentfix', overwrite=True, checksum=True)


def write_array_as_fits(outfile, data, header=None):
    """
    Helper function to take a 2D array object and save the data as a fits file.
    - If no header is supplied, then the basic info is written by Astropy.

    - This function circumvents sunpy.io in order to gain flexibility in the output
      data type and how the fits library is called.

    - I use it to write files of desired type for the Decurlog GPU deconvolution code.
    """
    # start a new fits file object
    hdu = astropy.io.fits.PrimaryHDU(data, header=header)

    # build the hdu list
    hdulist = astropy.io.fits.HDUList([hdu])

    # write the file
    hdulist.close()
    hdu.writeto(outfile, output_verify='silentfix', overwrite=True, checksum=True)


def write_array_as_compressed_fits(outfile, data, header=None, quantize_level=16.):
    """
    Helper function to take a 2D array object and save the data as a fits file.
    - If no header is supplied, then the basic info is written by Astropy.

    - This function circumvents sunpy.io in order to gain flexibility in the output
      data type and how the fits library is called.

    - I use it to write files of desired type for the Decurlog GPU deconvolution code.
    """
    # start a new fits file object
    hdu = astropy.io.fits.PrimaryHDU(data, header=header)

    # build the hdu list
    hdulist = astropy.io.fits.HDUList([hdu])

    # write out the file
    comp_hdu = astropy.io.fits.CompImageHDU(hdulist[0].data, header, quantize_level=quantize_level)
    hdulist.close()
    comp_hdu.writeto(outfile, output_verify='silentfix', overwrite=True, checksum=True)


def carrington_rotation_number_relative(time, lon):
    """
    A function that returns the decimal carrington rotation number for a spacecraft position
    that may not be at the same place at earth. In this case you know the carrington longitude
    of the spacecraft, and want to convert that to a decimal carrington number that is within
    +0.5 and -0.5 of the decimal rotation for the earth-based longitude.

    :param time: an astropy Time object indicating the time the position is known.
    :param lon: the carrington longitude of the spacecraft position.
    :return: the decimal_carrington number.
    """
    # get the decimal carrington number for Earth at this time
    cr_earth = sunpy.coordinates.sun.carrington_rotation_number(time)

    # convert that to the earth longitude (this should match sunpy.coordinates.sun.L0(time))
    cr0 = np.floor(cr_earth)
    lon_earth = np.mod((1 - (cr_earth - cr0)*360), 360)

    # compute the angular difference and the modulus
    diff = lon_earth - lon
    mod = np.mod(diff, 360.)

    # compute the fractional rotation offset, which depends on where the periodic boundary is.
    offset = 0.0
    if lon_earth < 180 and mod < 180 and diff < 0:
        offset = +1.0
    if lon_earth >= 180 and mod >= 180 and diff >= 0:
        offset = -1.0
    cr_now = cr0 + np.mod(1.0 - lon/360., 360.) + offset

    debug = False
    if debug:
        print('{: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f}'.format(lon, diff, mod, cr_now, cr_earth,
                                                                             cr_now - cr_earth))
        print(cr_earth, cr0, lon_earth, sunpy.coordinates.sun.L0(time).value, lon, cr_now)

    return cr_now


def construct_map_path_and_fname(base_dir, dtime, map_id, map_type, extension, inst=None, mkdir=True):
    """
    Wrapper to adapt construct_path_and_fname() for map files.
    - it returns the subdirectory path and filename
    """

    prefix = map_type
    postfix = 'MID' + str(map_id)
    maptype_base_dir = os.path.join(base_dir, map_type)
    # make the directory if needed
    if mkdir:
        # first check if the main directory exists
        if not os.path.isdir(base_dir):
            raise Exception('Base path does not exist! ' + base_dir)
            return None, None
        # check if the subdirectory exists
        if not os.path.isdir(maptype_base_dir):
            os.makedirs(maptype_base_dir)
    sub_dir, fname = construct_path_and_fname(maptype_base_dir, dtime, prefix, postfix, extension, inst=inst,
                                              mkdir=mkdir)

    return sub_dir, fname


def print_full_dataframe(df):
    """
    Helper function to print a pandas dataframe with NO truncation of rows/columns
    This will reset the defaults after, so it can be useful for inspection without annoying
    side effects in a script (pulled from a stack overflow example).
    """
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    pd.set_option('display.max_colwidth', None)
    print(df)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.max_colwidth')


def read_db_dir(dir_path):
    """
    Recursively walk through the directory and collect file dates and paths.

    Parameters
    ----------
    dir_path - character string
               The directory to be mapped.

    Returns
    -------
    A dataframe with date and relative-file path.
    """

    # intitialize empty dataframe
    out_df = pd.DataFrame(columns=['date', 'rel_path'])
    out_df['date'] = pd.to_datetime(out_df.date, utc=True)

    # create directory structure generator
    all_paths = os.walk(dir_path)

    for root, dirs, files in all_paths:
        if len(files) == 0:
            continue
        # check that this is not an 'index_files' directory
        if os.path.basename(root) == "index_files":
            # if so, these are not map/image files. Skip
            continue
        # first check for and remove '.DS_Store'
        valid_files = [x for x in files if (x != ".DS_Store") and (x != "README")]
        if len(valid_files) > 0:
            # loop through files and add to output data frame
            for filename in valid_files:
                # extract date (regex '\d{8}T\d{6}')
                date_str = re.search("\d{8}T\d{6}", filename).group()
                # construct relative path
                full_path = os.path.join(root, filename)
                rel_path = os.path.relpath(full_path, dir_path)
                new_row = pd.DataFrame(dict(date=[pd.to_datetime(date_str, utc=True), ], rel_path=[rel_path, ]))
                # add new row to output dataframe
                # out_df = out_df.append(new_row)
                out_df = pd.concat([out_df, new_row], axis=0, ignore_index=True)
        else:
            continue

    # sort dataframe by date (and reset row index)
    out_df.sort_values(by='date', inplace=True, ignore_index=True)

    return out_df


def gen_hipft_index(dir_path):
    available_maps = read_db_dir(dir_path)

    # save summary dataframe to file
    write_df = available_maps.rename(columns=dict(date='obs_datetime_utc', rel_path='map_path'))
    write_df.loc[:, 'target_datetime_utc'] = write_df.obs_datetime_utc.dt.round('H')
    # re-order columns and reset index
    write_df = write_df.loc[:, ['target_datetime_utc', 'obs_datetime_utc', 'map_path']]
    write_df.reset_index(drop=True, inplace=True)

    # add julian-days
    obs_astro_time = astro_time.Time(write_df.obs_datetime_utc)
    obs_jdays = obs_astro_time.jd
    # add new columns to dataframe
    jd_time_df = pd.DataFrame(dict(obs_jd=obs_jdays))
    write_df = pd.concat([jd_time_df, write_df], axis=1)
    # re-order columns
    write_df = write_df.loc[:, ['target_datetime_utc', 'obs_datetime_utc', 'obs_jd', 'map_path']]

    return write_df


def gen_hipft_index_old(dir_path):
    available_maps = read_db_dir(dir_path)

    # update column headers and add linux days
    write_df = available_maps.rename(columns=dict(date='hmi_datetime', rel_path='map_path'))
    write_df.loc[:, 'target_datetime'] = write_df.hmi_datetime.dt.round('H')
    # re-order columns and reset index
    write_df = write_df.loc[:, ['target_datetime', 'hmi_datetime', 'map_path']]
    write_df.reset_index(drop=True, inplace=True)

    # add fractional days since unix-epoch
    target_datetime = write_df.target_datetime.dt.to_pydatetime()
    target_unix_seconds = [float(target_datetime[ii].strftime("%s")) for ii in range(len(target_datetime))]
    target_unix_days = [x / (60 * 60 * 24) for x in target_unix_seconds]
    hmi_datetime = write_df.hmi_datetime.dt.to_pydatetime()
    hmi_unix_seconds = [float(hmi_datetime[ii].strftime("%s")) for ii in range(len(hmi_datetime))]
    hmi_unix_days = [x / (60 * 60 * 24) for x in hmi_unix_seconds]
    # add new columns to dataframe
    unix_time_df = pd.DataFrame(dict(target_unix_days=target_unix_days, hmi_unix_days=hmi_unix_days))
    write_df = pd.concat([unix_time_df, write_df], axis=1)

    return write_df
