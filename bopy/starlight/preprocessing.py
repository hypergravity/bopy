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


old comments
------------
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
# ####### the usage of np.seterr #######
# A=np.array([0.0])
# C=np.array([1.0])
# print A.dtype
# with np.errstate(divide='ignore'):
#    B=C/A
# print B
# ######################################


def preprocess_ccm(wave, flux, flux_err, mask, Av, z,
                   wave_interp,
                   Rv=3.1, strict_ccm=False):
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
        wave, flux, flux_err, Av=Av, Rv=Rv, strict_ccm=strict_ccm, z=None)
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
    P_flux = PchipInterpolator(wave_dez, flux_dez, extrapolate=False)
    flux_interp = P_flux(wave_interp)
    flux_interp = np.where(np.isnan(flux_interp), 0., flux_interp)

    # flux_err
    P_flux_err = PchipInterpolator(wave_dez, flux_err_dez, extrapolate=False)
    flux_err_interp = P_flux_err(wave_interp)
    flux_err_interp = np.where(np.isnan(flux_err_interp), 0., flux_err_interp)

    # mask
    mask_interper = interp1d(wave_dez, mask.astype(np.int64), kind='nearest',
                             fill_value=99, bounds_error=False)
    mask_interp = mask_interper(wave_interp)

    # strange mask spectra
    if np.sum(mask_interp == 0) == 0:
        # this is a bad spectrum
        mask_median = np.median(mask_interp)
        ind2 = np.logical_or(mask_interp / mask_median > 1000,
                             mask_interp / mask_median < 0)
        mask_interp[ind2] = 99
        mask_interp[np.logical_not(ind2)] = 0
    else:
        # this is a good spectrum
        mask_interp[mask_interp > 2] = 2

    return wave_interp, flux_interp, flux_err_interp, mask_interp


def preprocess_write_spectrum(wave, flux, flux_err, mask, outpath,
                              header=False, verbose=False, **kwargs):
    """ write spectra into STARLIGHT form

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
    # construct data table
    col_wave = Column(data=wave, name='wave', dtype=np.float64)
    col_flux = Column(data=flux, name='flux', dtype=np.float64)
    col_flux_err = Column(data=flux_err, name='flux_err', dtype=np.float64)
    col_mask = Column(data=mask, name='flag', dtype=np.int32)
    spec_reduc = Table(data=[col_wave, col_flux, col_flux_err, col_mask])

    # check basedir existence
    outpath_dirname = os.path.dirname(outpath)
    outpath_dirname2 = os.path.dirname(outpath_dirname)
    if not os.path.exists(outpath_dirname):
        # if this dir not existed
        try:
            assert os.path.exists(outpath_dirname2)
        except AssertionError:
            raise(AssertionError('@Cham: please create [%s] yourself!'
                                 % outpath_dirname2))
        os.mkdir(outpath_dirname)

    # write table
    if header:
        # ascii file with header
        spec_reduc.write(outpath, format='ascii')
    else:
        # ascii file with no header (default)
        spec_reduc.write(outpath, format='ascii.no_header')

    if verbose:
        print('@Cham: [write STARLIGHT input spectrum] %s ' % outpath)


def preprocess_sdss_dr10(Rv=3.1, strict_ccm=False,
                         filesource='sdss_dr10',
                         plate=2238, mjd=52525, fiberid=1,
                         Av=1., z=0.01,
                         wave_interp=None,
                         filepath_input='',
                         filepath_output='',
                         verbose=False):
    """ pre-processing a single SDSS DR10 spectrum

    Parameters
    ----------
    Rv: float
        the Rv parameter in CCM extinction law

    filesource: string
        the file source of spectra

    plate, mjd, fiberid: int
        the ID of SDSS spectra

    Av: float
        the magnitude of fore-ground extinction

    z: float
        the redshift of the galaxy

    wave_interp: array-like
        the wavelength array which interpolated into

    filepath_input: string
        the filepath of raw SDSS spectrum

    filepath_outpath: string
        the filepath of output pre-precessed spectra

    """

    # read spectrum
    spec = read_spectrum(filepath_input, filesource=filesource)

    # preprocess: de-red & de-z & re-sample
    wave, flux, flux_err, mask = preprocess_ccm(
        spec['wave'].data,
        spec['flux'].data,
        spec['flux_err'].data,
        spec['or_mask'].data,
        Av=Av,
        z=z,
        wave_interp=wave_interp,
        Rv=Rv,
        strict_ccm=strict_ccm)

    # preprocess: write spectrum
    preprocess_write_spectrum(wave, flux, flux_err, mask,
                              filepath_output,
                              verbose=verbose)


def preprocess_sdss_dr10_script(n_test, parallel=False):
    """ This is a script for pre-process spectra for STARLIGHT input

    Parameters
    ----------
    n_test: int
        number of spectra that will be run for test

    parallel: bool
        if True, use joblib to parallelize the code

    Example
    -------
    >>> from bopy.starlight.preprocessing import preprocess_sdss_dr10_script
    >>> ??preprocess_sdss_dr10_script()
    >>> # then copy them to a script & change parameters/settings

    """
    from bopy.starlight.preprocessing import preprocess_sdss_dr10
    from joblib import Parallel, delayed

    # settings
    Rv = 3.1
    strict_ccm = True

    # load catalog
    catalogpath = r'/pool/sdss/DR10/catalog/sdss.csv'
    catalog = Table.read(catalogpath)

    plate = catalog['plate'].data
    mjd = catalog['mjd'].data
    fiberid = catalog['fiberID'].data
    Av = catalog['ebv'].data / Rv
    z = catalog['z'].data

    # filepath for every spectra
    dir_input = r'/pool/sdss/DR10/spec'
    dir_output = r'/pool/sdss/DR10/spec_starlight'
    filepath_input = sdss_filepath(
        plate, mjd, fiberid, dir_input, extname='.fits')
    filepath_output = sdss_filepath(
        plate, mjd, fiberid, dir_output, extname='.dat')

    filesource = 'sdss_dr10'

    # wavelength array for interpolation
    wave_interp = np.arange(3400., 7100., 1.)

    # loop for every spectra
    if not parallel:
        # sequential version
        for i in xrange(n_test):
            preprocess_sdss_dr10(Rv=3.1, strict_ccm=strict_ccm,
                                 filesource=filesource,
                                 plate=plate[i], mjd=mjd[i], fiberid=fiberid[i],
                                 Av=Av[i], z=z[i],
                                 wave_interp=wave_interp,
                                 filepath_input=filepath_input[i],
                                 filepath_output=filepath_output[i],
                                 verbose=True)
    else:
        # parallel version
        Parallel(n_jobs=24, verbose=10)(delayed(preprocess_sdss_dr10)
            (Rv=3.1, strict_ccm=strict_ccm,
             filesource=filesource,
             plate=plate[i], mjd=mjd[i], fiberid=fiberid[i],
             Av=Av[i], z=z[i],
             wave_interp=wave_interp,
             filepath_input=filepath_input[i],
             filepath_output=filepath_output[i],
             verbose=True)
                                         for i in xrange(n_test))
    return


def preprocess_sdss_dr10_script_print():
    s = """
    from bopy.starlight.preprocessing import preprocess_sdss_dr10
    from joblib import Parallel, delayed

    # settings
    Rv = 3.1
    strict_ccm = True

    # load catalog
    catalogpath = r'/pool/sdss/DR10/catalog/sdss.csv'
    catalog = Table.read(catalogpath)

    plate = catalog['plate'].data
    mjd = catalog['mjd'].data
    fiberid = catalog['fiberID'].data
    Av = catalog['ebv'].data / Rv
    z = catalog['z'].data

    # filepath for every spectra
    dir_input = r'/pool/sdss/DR10/spec'
    dir_output = r'/pool/sdss/DR10/spec_starlight'
    filepath_input = sdss_filepath(
        plate, mjd, fiberid, dir_input, extname='.fits')
    filepath_output = sdss_filepath(
        plate, mjd, fiberid, dir_output, extname='.dat')

    filesource = 'sdss_dr10'

    # wavelength array for interpolation
    wave_interp = np.arange(3400., 7100., 1.)

    # loop for every spectra
    if not parallel:
        # sequential version
        for i in xrange(n_test):
            preprocess_sdss_dr10(Rv=3.1, strict_ccm=strict_ccm,
                                 filesource=filesource,
                                 plate=plate[i], mjd=mjd[i], fiberid=fiberid[i],
                                 Av=Av[i], z=z[i],
                                 wave_interp=wave_interp,
                                 filepath_input=filepath_input[i],
                                 filepath_output=filepath_output[i],
                                 verbose=True)
    else:
        # parallel version
        Parallel(n_jobs=24, verbose=10)(delayed(preprocess_sdss_dr10)
            (Rv=3.1, strict_ccm=strict_ccm,
             filesource=filesource,
             plate=plate[i], mjd=mjd[i], fiberid=fiberid[i],
             Av=Av[i], z=z[i],
             wave_interp=wave_interp,
             filepath_input=filepath_input[i],
             filepath_output=filepath_output[i],
             verbose=True)
                                         for i in xrange(n_test))
    return
    """
    print s
    return


if __name__ == '__main__':
    # preprocess_sdss_dr10_script(1000, True)
    preprocess_sdss_dr10_script_print()

