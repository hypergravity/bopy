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
            wave_start, wave_stop, R=R, over_sample=over_sample)
    elif isfunction(R):
        # if R is a function
        return _generate_wave_array_R_func(
            wave_start, wave_stop, R=R, over_sample=over_sample, wave_test_step=wave_test_step)


# ########################################################################### #
# ################  generate wave array given delta_lambda  ################# #
# ########################################################################### #


def generate_wave_array_delta_lambda_fixed(wave_start, wave_stop,
                                           delta_lambda,
                                           over_sample=1.):
    """ generate a wavelength array matching the given delta_lambda (fixed) """
    return np.arange(wave_start, wave_stop, delta_lambda / over_sample)


def generate_wave_array_delta_lambda_func(wave_start, wave_stop,
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


# ########################################################################### #
# #################  spectral R (FWHM) and R_max (FWHM_min) ################# #
# ########################################################################### #


def find_R_for_wave_array(wave):
    """ find the R of wavelength array (sampling resolution array) """
    wave = wave.flatten()
    wave_diff = np.diff(wave)
    wave_diff_ = (wave_diff[1:] + wave_diff[:-1]) / 2.
    R = np.hstack((
        wave[0]/wave_diff[0], wave[1:-1]/wave_diff_, wave[-1]/wave_diff[-1]))
    return R


def find_R_max_for_wave_array(wave):
    """ find the maximum sampling resolution of a given wavelength array """
    return np.max(find_R_for_wave_array(wave))


def find_appropriate_R_interp_for_convolution(R_hi, R_hi_spec):

    R_hi_spec_R_hi = R_hi_spec / R_hi

    if R_hi_spec_R_hi > OVER_SAMPLING:
        # already over sampled
        R_interp = R_hi_spec
    elif R_hi_spec_R_hi >= 1.:
        # already over sampled
        R_interp = OVER_SAMPLING * R_hi_spec
    else:
        R_interp = OVER_SAMPLING * R_hi

    return R_interp


def find_gaussian_kernel_fwhm(R_hi, R_lo, return_type='fwhm'):
    assert R_hi > R_lo
    fwhm_hi = resolution2fwhm(R_hi)
    fwhm_lo = resolution2fwhm(R_lo)
    fwhm = np.sqrt(fwhm_lo**2. - fwhm_hi**2.)
    if return_type == 'fwhm':
        return fwhm
    elif return_type == 'R':
        return fwhm2resolution(fwhm)






def fwhm2sigma(fwhm):
    return fwhm / (2. * np.sqrt(2. * np.log(2.)))


def sigma2fwhm(sigma):
    return 2. * np.sqrt(2. * np.log(2.)) * sigma


def normalized_gaussian_array(x, b=0., c=1.):
    # a = 1. / (np.sqrt(2*np.pi) * c)
    ngs_arr = np.exp(-(x-b)**2./(2.*c**2.))
    return ngs_arr / np.sum(ngs_arr)


def generate_gaussian_kernel_array(fwhm_pixel_num, array_length):
    sigma_pixel_num = fwhm2sigma(fwhm_pixel_num)
    array_length = np.fix(array_length)
    if array_length % 2 == 0:
        array_length += 1
    xgs = np.arange(array_length) - (array_length-1)/2.
    return normalized_gaussian_array(xgs, b=0., c=sigma_pixel_num)


def conv_spec(spec, R_hi, R_lo, R_interp=None, wave_new=None,
              wave_new_oversample=5., verbose=True):
    """ to convolve high-R spectrum to low-R spectrum

    Parameters
    ----------

    spec: spec object (Table with 'wave' and 'flux' column)

    R_hi: high R
        original resolution

    R_lo: low R
        target resolution

    R_interp: float
        interpolation resolution under which convolution takes place

    wave_new: float or vector
        if float: this specifies the over-sample rate
        if voctor: this specifies the new wave array

    verbose: bool
        if True, print the details on the screen

    Returns
    -------
    spec: bopy.spec.spec.Spec
        Spec, based on astropy.table.Table class

    """

    if verbose:
        start = datetime.datetime.now()
    wave_max = np.max(spec['wave'])
    wave_min = np.min(spec['wave'])

    # R_hi, R_lo, R_interp, Resamp_hi, Resamp_lo
    # 1. find R_interp
    if verbose:
        print '--------------------------------------------------------------'
        print '@Cham: determining R_interp ...'
    if R_interp is None:
        # need to find an R_interp
        R_hi_specmax = find_spec_max_R(spec['wave'])
        R_interp = find_appropriate_R_interp_for_convolution(R_hi, R_hi_specmax)
    assert R_interp >= R_hi

    # 2. interpolate high resolution spectra (over-resample)
    if verbose:
        print '@Cham: interpolating orignal spectrum to R_interp ...'
    wave_interp = generate_wave_array_R(wave_min, wave_max, R_interp)
    flux_interp = pchip_interpolate(spec['wave'], spec['flux'], wave_interp)

    # under this R_interp, what kind of gaussian kernel do I need?
    # 3. calculate gaussian kernel array (under this R_interp, R_gaussian_kernel)
    if verbose:
        print '@Cham: generating gaussian kernel array ...'
    R_gaussian_kernel = find_gaussian_kernel_fwhm(R_hi, R_lo, 'R')
    fwhm_pixel_num = R_interp / R_gaussian_kernel
    gs_array = generate_gaussian_kernel_array(fwhm_pixel_num, KERNEL_LENGTH_FWHM*fwhm_pixel_num)
    gs_array_len = len(gs_array)
    gs_array_arm = (gs_array_len-1)/2.

    # 4. convolution
    if verbose:
        print '@Cham: spectrum convolution start ...'
        print '@Cham: R_hi:              %.2f' % R_hi
        print '@Cham: R_hi_specmax:      %.2f' % R_hi_specmax
        print '@Cham: R_interp:          %.2f' % R_interp
        print '@Cham: R_lo:              %.2f' % R_lo
        print '@Cham: R_gaussian_kernel: %.2f' % R_gaussian_kernel
        print '@Cham: convolution array length: [%d, %d]'\
              % (len(flux_interp), gs_array_len)
        print '@Cham: target spectrum R_lo: %.2f' % R_lo
        print '@Cham: spectrum convolving ...'
        print '@Cham: estimated convolution time: %.2f seconds ...' \
              % (len(flux_interp)*gs_array_len/55408./657.*0.05)
        start_conv = datetime.datetime.now()
    convolved_flux = np.convolve(flux_interp, gs_array)[gs_array_arm:-gs_array_arm]
    if verbose:
        stop_conv = datetime.datetime.now()
        print '@Cham: spectrum convolution complete ...'
        print '@Cham: convolution time spent: %.2f seconds' \
              % (stop_conv-start_conv).total_seconds()

    # 5. find new wave array
    if wave_new is None:
        # need to find new wave
        # default: 5 times over-sample
        if verbose:
            print '@Cham: using default 5 times over-sample wave array ...'
        wave_new = generate_wave_array_R(wave_interp[0], wave_interp[-1], wave_new_oversample*R_lo)
    elif np.isscalar(wave_new):
        # wave_new specifies the new wave array resampling rate
        # default is 5. times over-sample
        if verbose:
            print '@Cham: using user-specified %.2f times over-sample wave array ...' % wave_new
        wave_new = generate_wave_array_R(wave_interp[0], wave_interp[-1], wave_new*R_lo)
    else:
        if verbose:
            print '@Cham: using user-specified wave array ...'

    # 6. interpolate convolved flux to new wave array
    if verbose:
        print '@Cham: interpolating convolved spectrum to new wave array ...'
    flux_new = pchip_interpolate(wave_interp, convolved_flux, wave_new)
    if verbose:
        stop = datetime.datetime.now()
        print '@Cham: total time spent: %.2f seconds' % (stop-start).total_seconds()
        print '--------------------------------------------------------------'

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
    spec_ = conv_spec(spec, 2000, 500, verbose=True)
    spec_.pprint()

    # 3.plot results
    fig = plt.figure()
    plt.plot(spec['wave'], spec['flux'])
    plt.plot(spec_['wave'], spec_['flux'], 'r')
    fig.show()
    print '@Cham: test OK ...'


if __name__ == '__main__':
    test_bc03_degrade_to_R500()