# -*- coding: utf-8 -*-
"""
@author: cham
Created on Fri Aug 14 15:24:20 2015
"""

# import aplpy
# from astropy.table import Table
# from astropy.coordinates import Galactic, SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from reproject import reproject_from_healpix, reproject_interp, reproject_to_healpix
import matplotlib.pyplot as plt
from healpy import write_map


def image_reproject_from_healpix_to_file(source_image_hdu, target_image_hdu_header, filepath=None):
    """ reproject from healpix image to normal wcs image

    :param source_image_hdu:        the HDU object of source image (healpix)
    :param target_image_hdu_header: the HDU object of target image (wcs)
    :param filepath:                the output file path
    :return:                        array, footprint
    """
    array, footprint = reproject_from_healpix(source_image_hdu, target_image_hdu_header)
    if filepath is not None:
        # write file
        fits.writeto(filepath, array, target_image_hdu_header, clobber=True)  # clobber=OVERWRITE
    else:
        # return array & footprint
        return array, footprint


def image_reproject_to_healpix_to_file(array, target_image_hdu_header, coordsys='galactic', filepath=None):
    """reproject image array to healpix image and write file

    :param array:                   image data
    :param target_image_hdu_header: the HDU object of
    :param coordsys:                target coordinate system --> {'galactic', 'equatorial'}
    :param filepath:                the output file path
    :return:                        array, footprint (only if filepath=None)
    """
    array, footprint = reproject_to_healpix((array, target_image_hdu_header), coordsys)
    if filepath is not None:
        # write file
        write_map(filepath, array)
    else:
        # return array & footprint
        return array, footprint


def image_reproject_wcs_to_file(source_image_hdu, target_image_hdu_header, filepath=None):
    """reproject one wcs image to the wcs of another image

    :param source_image_hdu:        the HDU object of source image
    :param target_image_hdu_header: the header object of target image
    :param filepath:                the output file path
    :return:
    """
    array, footprint = reproject_interp(source_image_hdu, target_image_hdu_header)
    if filepath is not None:
        # write file
        fits.writeto(filepath, array, target_image_hdu_header, clobber=True)  # clobber=OVERWRITE
    else:
        # return array & footprint
        return array, footprint


def image_query(data_type, target_image_hdu_header, filepath=None):
    """query image from stored data

    :param data_type:               data type of the queried source image
    :param target_image_hdu_header: target header
    :param filepath:                output file path
    :return:                        array, footprint
    """
    data_type_coll = {
        'haslam_408':   ('healpix', 1, '/pool/maps/LAMBDA/haslam408/haslam408_dsds_Remazeilles2014_ns2048.fits'),
        'planck_30':    ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/LFI_SkyMap_030_1024_R2.01_full.fits'),
        'planck_44':    ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/LFI_SkyMap_044_1024_R2.01_full.fits'),
        'planck_70':    ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/LFI_SkyMap_070_2048_R2.01_full.fits'),
        'planck_100':   ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/HFI_SkyMap_100_2048_R2.00_full.fits'),
        'planck_143':   ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/HFI_SkyMap_143_2048_R2.00_full.fits'),
        'planck_217':   ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/HFI_SkyMap_217_2048_R2.00_full.fits'),
        'planck_353':   ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/HFI_SkyMap_353_2048_R2.00_full.fits'),
        'planck_545':   ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/HFI_SkyMap_545_2048_R2.00_full.fits'),
        'planck_857':   ('healpix', 1, '/pool/skysurvey/PLANCK/PR2/HFI_SkyMap_857_2048_R2.00_full.fits'),
        'planck_dust':  ('healpix', 1,
                         '/pool/skysurvey/PLANCK/PR2/com/COM_CompMap_ThermalDust-commander_2048_R2.00.fits')}
    assert data_type in data_type_coll.keys()
    data_type = 'planck_dust'
    map_type, map_hdu, map_path = data_type_coll[data_type]
    if map_type == 'healpix':
        if filepath is not None:
            image_reproject_from_healpix_to_file(fits.open(map_path)[map_hdu],
                                                 target_image_hdu_header, filepath=filepath)
        else:
            array, footprint = image_reproject_from_healpix_to_file(fits.open(map_path)[map_hdu],
                                                                    target_image_hdu_header, filepath=None)
            return array, footprint
    else:
        return None, None


# -------------------------------------------------------------------------
if __name__ == '__main__':
    target_header = fits.open('/pool/mmt/2015b/wise/ngc_663/w1/mosaic_bm/mosaic.fits')[0].header
    target_wcs = WCS(target_header)
    # haslam408 = fits.open('/pool/maps/LAMBDA/haslam408/haslam408_dsds_Remazeilles2014_ns2048.fits')[1]
    # haslam408 = fits.open('/pool/maps/LAMBDA/IRIS/IRIS_nohole_1_2048.fits')[1]
    # hdu_in = fits.open('/pool/MMT/2015b/iras/b1/mosaic16/mosaic.fits')[0]
    # array, footprint = image_reproject_wcs_to_file(hdu_in, target_header)
    array_, footprint_ = image_query('planck_857', target_header)
    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1, projection=target_wcs)
    ax1.imshow(array_, origin='lower', vmin=1, vmax=5000)
    ax1.coords.grid(color='white')
    ax1.coords['ra'].set_axislabel('Right Ascension')
    ax1.coords['dec'].set_axislabel('Declination')
    fig.canvas.draw()
