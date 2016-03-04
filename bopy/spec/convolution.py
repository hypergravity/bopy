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

import numpy as np
import datetime
import astropy.constants as const
import astropy.units as u
from scipy.interpolate import pchip_interpolate
from bopy.spec.spec import spec_quick_init

OVER_SAMPLING = 10.
# threshold for R_hi_spec/R_hi
# if R_hi_spec > OVER_SAMPLING*R_hi, then R_interp = R_hi_spec
# else R_interp = OVER_SAMPLING * np.max([R_hi_spec, R_hi])

KERNEL_LENGTH_FWHM = 4.24
# kernel array length is 4.24 times FWHM, i.e., 10 sigma


def resolution2fwhm(R):
    assert R > 0.
    return const.c.value/R


def fwhm2resolution(fwhm):
    assert fwhm > 0.
    return const.c.value/fwhm


def find_gaussian_kernel_fwhm(R_hi, R_lo, return_type='fwhm'):
    assert R_hi > R_lo
    fwhm_hi = resolution2fwhm(R_hi)
    fwhm_lo = resolution2fwhm(R_lo)
    fwhm = np.sqrt(fwhm_lo**2. - fwhm_hi**2.)
    if return_type == 'fwhm':
        return fwhm
    elif return_type == 'R':
        return fwhm2resolution(fwhm)


def generate_wave_array_step(wave_start, wave_stop, wave_step):
    return np.arange(wave_start, wave_stop, wave_step)


def generate_wave_array_R(wave_start, wave_stop, r):
    r_ = r - .5
    # initial guess
    wave_step_min = wave_start / r_
    wave_step_guess = np.zeros((wave_stop-wave_start)/wave_step_min)
    wave_guess = np.zeros_like(wave_step_guess)
    wave_step_guess[0] = wave_step_min
    wave_guess[0] = wave_start
    # iterate for real
    for i in np.arange(len(wave_guess)-1)+1:
        wave_guess[i] = wave_guess[i-1] + wave_step_guess[i-1]
        wave_step_guess[i] = wave_guess[i] / r_
    return wave_guess[
        np.logical_and(wave_guess >= wave_start, wave_guess <= wave_stop)]


def find_spec_max_R(wave):
    wave_diff = np.diff(wave)
    return np.max(2. * wave[1:-1] / (wave_diff[1:] + wave_diff[:-1]))


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