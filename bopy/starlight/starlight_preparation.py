# -*- coding: utf-8 -*-
"""

Author
------
Bo Zhang

Email
-----
bozhang@nao.cas.cn

Created on
----------
- Fri Jun 15 15:37:00 2016

Modifications
-------------
- Fri Jun 15 15:37:00 2016    re-format code

Aims
----
- to pre-process spectra in order to meet the needs of STARLIGHT


How many places should I modify if I run this code on my computer?
0. catalogpath, catalog column names
1. reddening parameters: R_V, z, strict_ccm
2. working directory: inputdir, outputdir, datarelease
3. effective wavelength for interpolation: wavei
4. rules for generating filepath: filepath, outpath
of course you should change your n_jobs according to your computer configuration

"""

import os
import string
import glob
import sys
import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d

import astropy.io.ascii as ascii
from astropy.table import Table, Column, vstack
from astropy.io import fits

from joblib import Parallel, delayed

from .deredden import dered_ccm
from ..spec.read_spectrum import read_spectrum
from ..spec.lamost import sdss_filepath

np.seterr(divide='ignore', invalid='ignore')


def preprocess_ccm(wave, flux, flux_err, mask, Av, z,
                         wave_interp,
                         Rv=3.1, verbose=False):
    """ pre-processing of spectra to meet the needs of STARLIGHT

    Parameters
    ----------
    wave: numpy.ndarray
        wavelength array

    flux: numpy.ndarray
        flux array

    flux_err: numpy.ndarray
        flux error array

    mask: numpy.ndarray
        mask array

    Av: float
        the Galactic fore-ground extinction

    z: float
        the redshift of galaxy

    wave_interp: numpy.ndarray
        the wavelength array interpolated to

    verbose:

    Returns:

    """

    # 1. flux calibration
    flux *= 1E-17
    flux_err *= 1E-17

    # 2. de-reddening
    # Compute host reddening law for -effective observed wavelength wave,
    # at redshift z, assuming Rv.
    flux_dered, flux_err_dered = dered_ccm(
        wave, flux, flux_err, Av=Av, Rv=Rv, z=None, strict_ccm=False)
    # Al = entry['ebv'] * (a * R_V + b) --> WHY?

    # 3. de-redshift
    # wave corrected by a factor of 1/(1+z)
    wave_dez = wave / (1 + z)
    # flux corrected by a factor of (1+z)
    flux_dez = flux_dered * (1 + z)
    # flux corrected by a factor of (1+z)
    flux_err_dez = flux_err_dered * (1 + z)

    # 4. interpolate spectrum
    ###########################################################################
    # kind : str or int, optional
    # ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’
    # where ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline
    # interpolation of first, second or third order)
    # or as an integer specifying the order of the spline interpolator to use.
    # Default is ‘linear’.
    ###########################################################################

    # flux
    P_flux = PchipInterpolator(wave_dez, flux_dez, extrapolate=0)
    flux_interp = P_flux(wave_interp)

    # nois
    P_flux_err = PchipInterpolator(wave_dez, flux_err_dez, extrapolate=99.)
    flux_err_interp = P_flux_err(wave_interp)

    # mask
    mask_interp = interp1d(wave_dez, mask, kind='nearest')

    # strange mask spectra
    if np.sum(mask == 0) == 0:
        # this is a bad spectrum
        mask_median = np.median(mask)
        ind2 = np.logical_or(mask / mask_median > 1000, mask / mask_median < 0)
        mask[ind2] = 2
        mask[np.logical_not(ind2)] = 0
    else:
        # this is a good spectrum
        ind0 = mask == 0
        mask[np.logical_not(ind0)] = 2

    return wave_interp, flux_interp, flux_err_interp, mask_interp


def preprocess_write_spectrum(wave, flux, flux_err, mask, outpath,
                              header=False, verbose=False):
    """" write spectra into STARLIGHT form

    Parameters
    ----------
    wave: numpy.ndarray
        wavelength array

    flux: numpy.ndarray
        flux array

    flux_err: numpy.ndarray
        flux error array

    mask: numpy.ndarray
        mask array

    header: bool
        whether include headers in output spectrum

    verbose: bool
        whether print info during processing

    """
    spec_reduc = Table([wave, flux, flux_err, mask],
                       names=['wave', 'flux', 'flux_err', 'flag'])
    if header:
        # ascii file with header
        spec_reduc.write(outpath, format='ascii')
    else:
        # ascii file with no header (default)
        spec_reduc.write(outpath, format='ascii.no_header')

    if verbose:
        print('@Cham: [write STARLIGHT input spectrum] %s ' % outpath)



def preprocess():
    pass


def preprocess_sdss_dr10():
    Rv = 3.1
    strict_ccm = True

    # load catalog
    catalogpath = r'/pool/SDSS.DR10/catalog/sdss.csv'
    catalog = Table.read(catalogpath)

    Av = catalog['ebv'].data/Rv
    plate = catalog['plate'].data
    mjd = catalog['mjd'].data
    fiberid = catalog['fiberID'].data
    z = catalog['z']

    dir_input = r'/pool/SDSS/DR12/spectra'
    dir_output = r'/pool/SDSS/DR12/starlight'

    filepath_input = sdss_filepath(
        plate, mjd, fiberid, dir_input, extname='.fits')
    filepath_output = sdss_filepath(
        plate, mjd, fiberid, dir_input, extname='.dat')

    filesource = 'sdss_dr10'
    wave_interp = np.arange(3000., 8500., 1.)

    for i in xrange(len(z)):
        spec = read_spectrum(filepath_input[i], filesource=filesource)
        wave, flux, flux_err, mask = preprocess_ccm(
            spec['wave'].data,
            spec['flux'].data,
            spec['flux_err'].data,
            spec['or_mask'].data,
            Av[i],
            z=z[i],
            wave_interp=wave_interp,
            Rv=Rv,
            strict_ccm=strict_ccm,
            verbose=False)
        preprocess_write_spectrum(wave, flux, flux_err, mask,
                                  filepath_output[i], verbose=True)


# def starlight_preparation(ShowOnMonitor):
#
#     # joblib Parallel execution
#     Parallel(n_jobs=24, verbose=100)(delayed(preprocess()) \
#                                          (entry, ShowOnMonitor) for entry in
#                                      catalog)


if __name__ == '__main__':
    preprocess_sdss_dr10(ShowOnMonitor=False)

'''**************************************************************************
How many places should I modify if I run this code on my computer?
0. catalogpath, catalog column names
1. reddening parameters: R_V, z, strict_ccm
2. working directory: inputdir, outputdir, datarelease
3. effective wavelength for interpolation: wavei
4. rules for generating filepath: filepath, outpath
of course you should change your n_jobs according to your computer configuration
*****************************************************************************
'''

# A=np.array([0.0])
# C=np.array([1.0])
# print A.dtype
# with np.errstate(divide='ignore'):
#    B=C/A
# print B
