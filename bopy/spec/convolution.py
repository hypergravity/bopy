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
- Fri Feb  26 21:44:00 2016     convolution

Modifications
-------------
-

Aims
----
- degrade high resolution spectra to low resolution spectra

"""

# ###########################################################################
# Enlightening example:
#
# for one pixel:
# R_hi = 2000 @5000A    FWHM_hi =  2.5A
# R_lo =  500 @5000A    FWHM_lo = 10.0A
# FWHM_GK = sqrt(10.0**2 - 2.5**2) = 9.68A
#
# for next pixel:
# R_hi = 1000 @5000A    FWHM_hi =  5.0A
# R_lo =  500 @5000A    FWHM_lo = 10.0A
# FWHM_GK = sqrt(10.0**2 - 5.0**2) = 8.66A
#
# Therefore, to keep the number of pixels of Gaussian Kernels the same,
# we have to adjust the FWHM_interp (and thus the R_interp).
#
# For pixels 2000->500, the FWHM_GK = 9.68 (wider GK),
# while for pixels 1000-> 500, FWHM_GK = 8.66 (narrower GK),
# thus the R_interp should be 9.68/8.66=1.12 times higher,
# i.e., the delta_lambda should be  1./1.12 times smaller.
# ###########################################################################

import numpy as np
import datetime
import astropy.constants as const
import astropy.units as u
from inspect import isfunction
from scipy.interpolate import pchip_interpolate
from bopy.spec.spec import spec_quick_init

OVER_SAMPLING = 10.
# threshold for R_hi_spec/R_hi
# if R_hi_spec > OVER_SAMPLING*R_hi, then R_interp = R_hi_spec
# else R_interp = OVER_SAMPLING * np.max([R_hi_spec, R_hi])

KERNEL_LENGTH_FWHM = 4.24
# kernel array length is 4.24 times FWHM, i.e., 10 sigma


# ########################################################################### #
# ################  transformation between R and FWHM  ###################### #
# ########################################################################### #


def resolution2fwhm(R, wave=5000.):
    # assert R is positive
    if np.isscalar(R):
        assert R > 0.
    else:
        assert np.all(R > 0.)

    # assert wave is correct in shape
    assert np.isscalar(wave) or len(wave) == len(R)

    return wave / R


def fwhm2resolution(fwhm, wave=5000.):
    # assert fwhm is positive
    if np.isscalar(fwhm):
        assert fwhm > 0.
    else:
        assert np.all(fwhm > 0.)

    # assert wave is correct in shape
    assert np.isscalar(wave) or len(wave) == len(fwhm)

    return wave / fwhm


# ########################################################################### #
# ######################  generate wave array given R  ###################### #
# ########################################################################### #


def _generate_wave_array_R_fixed(wave_start, wave_stop, R=2000.,
                                over_sample=1.):
    """ generate wave array given a fixed R """
    R_ = over_sample * R - .5
    # determine wave_step_min
    wave_step_min = wave_start / R_
    # wave guess
    wave_step_guess = np.zeros((wave_stop-wave_start)/wave_step_min)
    wave_guess = np.zeros_like(wave_step_guess)
    wave_step_guess[0] = wave_step_min
    wave_guess[0] = wave_start
    # iterate for real
    for i in np.arange(1, len(wave_guess)):
        wave_guess[i] = wave_guess[i-1] + wave_step_guess[i-1]
        wave_step_guess[i] = wave_guess[i] / R_
    return wave_guess[
        np.logical_and(wave_guess >= wave_start, wave_guess <= wave_stop)]


def _generate_wave_array_R_func(wave_start, wave_stop, R=(lambda x: x),
                               over_sample=1., wave_test_step=1.):
    """ generate wave array given R as a function """
    # determine wave_step_min
    wave_test = np.arange(wave_start, wave_stop, wave_test_step)
    R_test = over_sample * R(wave_test)
    wave_step_min = np.min(wave_test / R_test)
    # wave guess
    wave_guess = np.zeros((wave_stop-wave_start)/wave_step_min)
    wave_guess[0] = wave_start
    # iterate for real # single side R !!!
    for i in np.arange(1, len(wave_guess)):
        wave_guess[i] = wave_guess[i-1] + \
                        wave_guess[i-1] / (over_sample * R(wave_guess[i-1]))
    return wave_guess[
        np.logical_and(wave_guess >= wave_start, wave_guess <= wave_stop)]


def generate_wave_array_R(wave_start, wave_stop, R=2000.,
                          over_sample=1., wave_test_step=1.):
    """ generate a wavelength array matching the given R

    Parameters
    ----------
    wave_start: float
        start from this wavelength
    wave_stop: float
        stop at this wavelength
    R: float or function
        specify a fixed R or specify R as a function of wavelength
    over_sample: float
        over-sampling rate, default is 1.
    wave_test_step:
        used to determing the wave_step_min

    Returns
    -------
    wave: array
        an array matching the given R

    Example
    -------
    >>> def R(x): return 0.2*x
    >>> wave_array_R = generate_wave_array_R(4000., 5000., R)

    """
    if np.isscalar(R):
        # if R is scalar
        return _generate_wave_array_R_fixed(
            wave_start, wave_stop, R=R,
            over_sample=over_sample)
    else:
        # if R is a function / Interpolator
        return _generate_wave_array_R_func(
            wave_start, wave_stop, R=R,
            over_sample=over_sample, wave_test_step=wave_test_step)


# ########################################################################### #
# ################  generate wave array given delta_lambda  ################# #
# ########################################################################### #


def _generate_wave_array_delta_lambda_fixed(wave_start, wave_stop,
                                            delta_lambda,
                                            over_sample=1.):
    """ generate a wavelength array matching the given delta_lambda (fixed) """
    return np.arange(wave_start, wave_stop, delta_lambda / over_sample)


def _generate_wave_array_delta_lambda_func(wave_start, wave_stop,
                                           delta_lambda=(lambda x: 1.),
                                           over_sample=1.,
                                           wave_test_step=1.):
    """ generate a wavelength array matching the given delta_lambda
        (specified as a function of wavelength) """
    wave_test = np.arange(wave_start, wave_stop, wave_test_step)
    delta_lambda_min = np.min(delta_lambda(wave_test))
    wave_guess = np.arange(wave_start, wave_stop, delta_lambda_min)
    for i in xrange(1, len(wave_guess)):
        wave_guess[i] = \
            wave_guess[i-1] + delta_lambda(wave_guess[i-1]) / over_sample
    return wave_guess[
        np.logical_and(wave_guess >= wave_start, wave_guess <= wave_stop)]


def generate_wave_array_delta_lambda(wave_start, wave_stop,
                                     delta_lambda=(lambda x: 1.),
                                     over_sample=1.,
                                     wave_test_step=1.):
    """ generate a wavelength array matching the given delta_lambda
        (delta_lambda given as a fixed number or a function)
    Parameters
    ----------
    wave_start: float
        where the wavelength starts
    wave_stop: float
        where the wavelength stops
    delta_lambda: float or function
        specifies the delta_lambda as a fixed number or a function of wavelength
    over_sample: float
        over-sampling
    wave_test_step: float
        tests for the smallest wave_guess step

    Returns
    -------
    wave_guess: array

    Example
    -------
    >>> def dl(x): return 0.002*x
    >>> wave_array_dl = generate_wave_array_delta_lambda(4000., 5000., dl)

    """
    if np.isscalar(delta_lambda):
        # if delta_lambda is scalar
        return _generate_wave_array_delta_lambda_fixed(
            wave_start, wave_stop, delta_lambda=delta_lambda,
            over_sample=over_sample)
    else:
        # if delta_lambda is a function / Interpolator
        return _generate_wave_array_delta_lambda_func(
            wave_start, wave_stop, delta_lambda=delta_lambda,
            over_sample=over_sample, wave_test_step=wave_test_step)


# ########################################################################### #
# ############## find spectral R (FWHM) and R_max (FWHM_min) ################ #
# ########################################################################### #


def find_R_for_wave_array(wave):
    """ find the R of wavelength array (sampling resolution array) """
    wave = wave.flatten()
    wave_diff = np.diff(wave)
    wave_diff_ = (wave_diff[1:] + wave_diff[:-1]) / 2.
    return np.hstack((
        wave[0]/wave_diff[0], wave[1:-1]/wave_diff_, wave[-1]/wave_diff[-1]))


def find_R_max_for_wave_array(wave):
    """ find the maximum sampling resolution of a given wavelength array """
    return np.max(find_R_for_wave_array(wave))


def find_delta_lambda_for_wave_array(wave):
    """ find the delta_lambda of wavelength array (delta_lambda array) """
    wave = wave.flatten()
    wave_diff = np.diff(wave)
    wave_diff_ = (wave_diff[1:] + wave_diff[:-1]) / 2.
    return np.hstack((wave_diff[0], wave_diff_, wave[-1]))


def find_delta_lambda_min_for_wave_array(wave):
    """ find the minimum delta_lambda of a given wavelength array """
    return np.min(find_delta_lambda_for_wave_array(wave))


# ########################################################################### #
# ###############################  find Rgk  ################################ #
# ########################################################################### #


def find_Rgk(R_hi=2000., R_lo=500., over_sample=1.):
    """ find Rgk as a function of wavelength

    Parameters
    ----------
    R_hi: float or funtion
        higher resolution (as a function of wavelength)
    R_lo: float or funtion
        lower resolution (as a function of wavelength)
    over_sample: float
        over-sampled resolution, default is 1.

    Returns
    -------
    Rgk: function
        Gaussian Kernel resolution as a function of wavelength

    """
    if np.isscalar(R_hi):
        # if R_hi is a fixed number
        R_hi_ = lambda x: R_hi
    else:
        R_hi_ = R_hi

    if np.isscalar(R_lo):
        # if R_lo is a fixed number
        R_lo_ = lambda x: R_lo
    else:
        R_lo_ = R_lo

    Rgk = lambda x: \
        (over_sample * x / np.sqrt((x/R_lo_(x))**2. - (x/R_hi_(x))**2.))
    return Rgk


# def find_appropriate_R_interp_for_convolution(R_hi, R_hi_spec):
#
#     R_hi_spec_R_hi = R_hi_spec / R_hi
#
#     if R_hi_spec_R_hi > OVER_SAMPLING:
#         # already over sampled
#         R_interp = R_hi_spec
#     elif R_hi_spec_R_hi >= 1.:
#         # already over sampled
#         R_interp = OVER_SAMPLING * R_hi_spec
#     else:
#         R_interp = OVER_SAMPLING * R_hi
#
#     return R_interp


# def find_gaussian_kernel_fwhm(R_hi, R_lo, return_type='fwhm'):
#     assert R_hi > R_lo
#     fwhm_hi = resolution2fwhm(R_hi)
#     fwhm_lo = resolution2fwhm(R_lo)
#     fwhm = np.sqrt(fwhm_lo**2. - fwhm_hi**2.)
#     if return_type == 'fwhm':
#         return fwhm
#     elif return_type == 'R':
#         return fwhm2resolution(fwhm)

# ########################################################################### #
# #########################  Gaussian Kernel ################################ #
# ########################################################################### #


def fwhm2sigma(fwhm):
    return fwhm / (2. * np.sqrt(2. * np.log(2.)))


def sigma2fwhm(sigma):
    return 2. * np.sqrt(2. * np.log(2.)) * sigma


def normalized_gaussian_array(x, b=0., c=1.):
    # a = 1. / (np.sqrt(2*np.pi) * c)
    ngs_arr = np.exp(- (x-b)**2. / (2.*c**2.))
    return ngs_arr / np.sum(ngs_arr)


def generate_gaussian_kernel_array(over_sample_Rgk, sigma_num):
    """ generate gaussian kernel array according to over_sample_Rgk

    Parameters
    ----------
    over_sample_Rgk: float
        over_sample rate
    sigma_num: float
        1 sigma of the Gaussian = sigma_num pixels

    Returns
    -------
    normalized gaussian array

    """
    sigma_pixel_num = fwhm2sigma(over_sample_Rgk)

    array_length = np.fix(sigma_num * sigma_pixel_num)
    if array_length % 2 == 0:
        array_length += 1
    array_length_half = (array_length-1) / 2.

    xgs = np.arange(- array_length_half, array_length_half + 1)
    return normalized_gaussian_array(xgs, b=0., c=sigma_pixel_num)


def conv_spec(spec,
              R_hi=2000.,
              R_lo=500.,
              over_sample_additional=3.,
              gaussian_kernel_sigma_num=8.,
              wave_new=None,
              wave_new_oversample=5.,
              verbose=True):
    """ to convolve high-R spectrum to low-R spectrum

    Parameters
    ----------
    spec: spec object (Table with 'wave' and 'flux' column)
        the original spectrum
    R_hi: float or function
        higher resolution
    R_lo: float or function
        lower R
    over_sample_additional: float
        additional over-sample rate
    gaussian_kernel_sigma_num: float
        the gaussian kernel width in terms of sigma
    wave_new: None or float or array
        if None: wave_new auto-generated using wave_new_oversample
        if float: this specifies the over-sample rate
        if voctor: this specifies the new wave_new array
    verbose: bool
        if True, print the details on the screen

    Returns
    -------
    spec: bopy.spec.spec.Spec
        Spec, based on astropy.table.Table class
    gk_len_half: float
        number of pixels in the head and tail which is bad

    """
    if verbose:
        start = datetime.datetime.now()
        print '---------------------------------------------------------------'
        print '@Cham: Welcome to the spectral convolution code developed by Bo Zhang (@NAOC) ...'

    # 1. re-format R_hi & R_lo
    assert R_hi is not None and R_lo is not None

    if np.isscalar(R_hi):
        R_hi_ = lambda x: R_hi
    else:
        R_hi_ = R_hi

    if np.isscalar(R_lo):
        R_lo_ = lambda x: R_lo
    else:
        R_lo_ = R_lo

    # 2. find Rgk
    Rgk = find_Rgk(R_hi_, R_lo_, over_sample=1.)

    # 3. find appropriate over_sample
    R_hi_specmax = find_R_max_for_wave_array(spec['wave'])
    R_hi_max = np.max(R_hi_(spec['wave']))
    over_sample = over_sample_additional * np.fix(np.max([
        R_hi_specmax/Rgk(spec['wave']), R_hi_max/Rgk(spec['wave'])]))

    # 4. find wave_interp & flux_interp
    if verbose:
        print '@Cham: interpolating orignal spectrum to wave_interp ...'
    wave_max = np.max(spec['wave'])
    wave_min = np.min(spec['wave'])
    wave_interp = generate_wave_array_R(wave_min, wave_max,
                                        Rgk, over_sample=over_sample)
    flux_interp = pchip_interpolate(spec['wave'], spec['flux'], wave_interp)

    # 5. generate Gaussian Kernel array
    if verbose:
        print '@Cham: generating gaussian kernel array ...'
    gk_array = generate_gaussian_kernel_array(over_sample,
                                              gaussian_kernel_sigma_num)
    gk_len = len(gk_array)
    gk_len_half = (gk_len - 1) / 2.

    # 6. convolution
    if verbose:
        print '@Cham: convolution ...'
        print '@Cham: estimated convolution time: %.2f seconds ...'\
              % (len(flux_interp)*gk_len/55408./657.*0.05)
    convolved_flux = np.convolve(flux_interp, gk_array)[gk_len_half:-gk_len_half]

    # 7. find new wave array
    if wave_new is None:
        # wave_new is None
        # default: 5 times over-sample
        if verbose:
            print '@Cham: using default 5 times over-sample wave array ...'
        wave_new = generate_wave_array_R(wave_interp[0], wave_interp[-1],
                                         R_lo, wave_new_oversample)
    elif np.isscalar(wave_new):
        # wave_new specifies the new wave array over_sampling_lo rate
        # default is 5. times over-sample
        if verbose:
            print '@Cham: using user-specified %.2f times over-sample wave array ...' % wave_new
        wave_new = generate_wave_array_R(wave_interp[0], wave_interp[-1],
                                         R_lo, wave_new)
    else:
        # wave_new specified
        if verbose:
            print '@Cham: using user-specified wave array ...'

    # 8. interpolate convolved flux to new wave array
    if verbose:
        print '@Cham: interpolating convolved spectrum to new wave array ...'
    flux_new = pchip_interpolate(wave_interp, convolved_flux, wave_new)
    if verbose:
        stop = datetime.datetime.now()
        print '@Cham: total time spent: %.2f seconds' % (stop-start).total_seconds()
        print '---------------------------------------------------------------'

    return spec_quick_init(wave_new, flux_new)


def test_bc03_degrade_to_R500():
    # test a BC03 population spectrum
    from astropy.table import Table
    import matplotlib.pyplot as plt
    from bopy.spec.spec import Spec
    # 1. read spectrum
    fp = '/home/cham/PycharmProjects/bopy/bopy/data/model_bc03/bc2003_hr_m42_chab_ssp_020.spec'
    data = np.loadtxt(fp)
    spec = Spec(data, names=['wave', 'flux'])
    print spec
    spec = spec.extract_chunk_wave_interval([[4000., 8000.]])[0]

    # 2.convolve spectum
    spec_ = conv_spec(spec, lambda x: 0.2*x, lambda x: 0.1*x,
                      over_sample_additional=3.,
                      gaussian_kernel_sigma_num=6.,
                      wave_new=None,
                      wave_new_oversample=5.,
                      verbose=False)
    spec_.pprint()

    print find_R_for_wave_array(spec_['wave'])
    # 3.plot results
    # fig = plt.figure()
    # plt.plot(spec['wave'], spec['flux'])
    # plt.plot(spec_['wave'], spec_['flux'], 'r')
    # fig.show()
    # fig.savefig(''')
    print '@Cham: test OK ...'


if __name__ == '__main__':
    test_bc03_degrade_to_R500()