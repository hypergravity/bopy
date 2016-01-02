# -*- coding: utf-8 -*-
"""

Author
------
- Bo Zhang

Email
-----
- bozhang@nao.cas.cn

Created on
----------
- Fri Jul  3 13:13:06 2015    read_spectrum

Modifications
-------------
- Wed Jul 29 21:46:00 2015    measure_line_index
- Fri Nov 20 10:16:59 2015    reformatting code

Aims
----
- to read LAMOST/SDSS spectra
- measure line index from spectra

"""

from astropy.io import fits
from astropy.table import Table, Column
import numpy as np
# from lmfit.models import LinearModel, GaussianModel
from lmfit.models import LinearModel, GaussianModel



def lamost_filepath(planid, mjd, spid, fiberid,
                    dirpath=r'/pool/lamost/dr2/spectra/fits/'):
    """generate file path of A LAMOST spectrum

    :param planid:  planid
    :param mjd:     mjd
    :param spid:    spid
    :param fiberid: fiberid
    :param dirpath: the root directory for storing spectra
    :return:        filepath
    """
    return ('%s%s/spec-%05d-%s_sp%02d-%03d.fits') % (dirpath, planid, mjd, planid, spid, fiberid)


def read_spectrum(filepath, filesource='auto'):
    """read SDSS/LAMOST spectrum

    Parameters
    ----------
    :param filepath:    input file path
    :param filesource:  'sdss_dr12'/'lamost_dr2'
    :return:
    """
    # auto-identify the spectrum origination
    if filesource == 'auto':
        telescope = fits.open(filepath)[0].header['TELESCOP']
        if telescope == 'SDSS 2.5-M':
            return read_spectrum(filepath, filesource='sdss_dr12')
        if telescope == 'LAMOST':
            return read_spectrum(filepath, filesource='lamost_dr2')
    # SDSS DR10/DR12 spectrum
    if filesource == 'sdss_dr10' or filesource == 'sdss_dr12':
        data = fits.open(filepath)
        specdata = Table(data[1].data)
        wave = Column(name='wave', data=np.power(10., specdata['loglam']))
        flux_err = Column(name='flux_err', data=(specdata['ivar']) ** -0.5)
        specdata.add_columns([wave, flux_err])
        return specdata
    # LAMOST DR2 spectrum
    if filesource == 'lamost_dr2':
        data = fits.open(filepath)
        specdata = Table(data[0].data.T)
        flux = Column(name='flux', data=specdata['col0'])
        ivar = Column(name='ivar', data=specdata['col1'])
        flux_err = Column(name='flux_err', data=(specdata['col1']) ** -0.5)
        wave = Column(name='wave', data=specdata['col2'])
        and_mask = Column(name='and_mask', data=specdata['col3'])
        or_mask = Column(name='or_mask', data=specdata['col4'])
        return Table([wave, flux, flux_err, ivar, and_mask, or_mask])
    return None


def measure_line_index(wave, flux, flux_err, mask,
                       line_info=None, num_refit=0, z=0.):
    """measure line index / line EW

    :param wave:        wavelength
    :param flux:        flux
    :param flux_err:    flux error
    :param mask:        andmask / ormask
    :param line_info:   information about spectral line (dict)
    :param num_refit:   number of refitting
    :param z:           redshift (only specify when z is large)
    :return:
    """
    try:
        ''' 0. import packages ------------------------------------------------
        '''
        from lmfit import minimize, Parameters, Model
        ''' 1. get line information -------------------------------------------
        '''
        # line_center = line_info['line_center']    # not used
        line_range = line_info['line_range']
        line_shoulder_left = line_info['line_shoulder_left']
        line_shoulder_right = line_info['line_shoulder_right']

        ''' 2. shift spectra to rest-frame ------------------------------------
        '''
        wave = np.array(wave)
        flux = np.array(flux)
        if z != 0:
            wave /= 1. + z

        ''' 3. estimate the local continuum -----------------------------------
        '''
        # shoulder wavelength range
        # -------------------------
        ind_shoulder = np.any([
            np.all([wave > line_shoulder_left[0],
                    wave < line_shoulder_left[1]], axis=0),
            np.all([wave > line_shoulder_right[0],
                    wave < line_shoulder_right[1]], axis=0)], axis=0)
        wave_shoulder = wave[ind_shoulder]
        flux_shoulder = flux[ind_shoulder]

        # integrated/fitted wavelength range
        # ----------------------------------
        ind_range = np.logical_and(wave > line_range[0], wave < line_range[1])
        wave_range = wave[ind_range]
        flux_range = flux[ind_range]
        # flux_err_range = flux_err[ind_range]  # not used
        mask_range = mask[ind_range]
        flux_err_shoulder = flux_err[ind_shoulder]
        # mask_shoulder = mask[ind_shoulder]    # not used

        ''' 4. linear model ---------------------------------------------------
        '''
        # from lmfit.models import LinearModel, GaussianModel
        mod_linear = LinearModel(prefix='mod_linear_')
        par_linear = mod_linear.guess(flux_shoulder, x=wave_shoulder)
        # to see the parameter names:
        # model_linear.param_names
        # {'linear_fun_intercept', 'linear_fun_slope'}
        out_linear = mod_linear.fit(flux_shoulder,
                                    par_linear,
                                    x=wave_shoulder,
                                    method='leastsq')
        ''' 4. estimate continuum ---------------------------------------------
        '''
        cont_shoulder = out_linear.best_fit
        noise_std = np.std(flux_shoulder / cont_shoulder)
        cont_range = mod_linear.eval(out_linear.params, x=wave_range)
        resi_range = 1 - flux_range / cont_range

        ''' 5.1 Integrated EW -------------------------------------------------
        '''
        # estimate EW_int
        # ---------------
        wave_diff = np.diff(wave_range)
        wave_step = np.mean(np.vstack([np.hstack([wave_diff[0], wave_diff]),
                                       np.hstack([wave_diff, wave_diff[-1]])]),
                            axis=0)
        EW_int = np.dot(resi_range, wave_step)

        # estimate EW_int_err
        # -------------------
        EW_int_err = np.std(np.dot(
            (resi_range.reshape(1, -1).repeat(100, axis=0) +
             np.random.randn(100, resi_range.size) * noise_std),
            wave_step))

        ''' 5.2 gaussian model ------------------------------------------------
        '''
        # estimate EW_fit
        # ---------------
        mod_gauss = GaussianModel(prefix='mod_gauss_')
        par_gauss = mod_gauss.guess(resi_range, x=wave_range)
        out_gauss = mod_gauss.fit(resi_range, par_gauss, x=wave_range)
        line_indx = {
            'SN_local_flux_err':        np.median(flux_shoulder / flux_err_shoulder),
            'SN_local_flux_std':        1. / noise_std,
            'flag_good_pixel':          np.all(mask_range == 0),
            'EW_int':                   EW_int,
            'EW_int_err':               EW_int_err,
            'mod_linear_slope':         out_linear.params[mod_linear.prefix + 'slope'].value,
            'mod_linear_slope_err':     out_linear.params[mod_linear.prefix + 'slope'].stderr,
            'mod_linear_intercept':     out_linear.params[mod_linear.prefix + 'intercept'].value,
            'mod_linear_intercept_err': out_linear.params[mod_linear.prefix + 'intercept'].stderr,
            'mod_gauss_amplitude':      out_gauss.params[mod_gauss.prefix + 'amplitude'].value,
            'mod_gauss_amplitude_err':  out_gauss.params[mod_gauss.prefix + 'amplitude'].stderr,
            'mod_gauss_center':         out_gauss.params[mod_gauss.prefix + 'center'].value,
            'mod_gauss_center_err':     out_gauss.params[mod_gauss.prefix + 'center'].stderr,
            'mod_gauss_sigma':          out_gauss.params[mod_gauss.prefix + 'sigma'].value,
            'mod_gauss_sigma_err':      out_gauss.params[mod_gauss.prefix + 'sigma'].stderr,
            'mod_gauss_amplitude_std':  np.nan,
            'mod_gauss_center_std':     np.nan,
            'mod_gauss_sigma_std':      np.nan}

        ''' 6. add noise and re-fit +++++++++++++++++++++++++++++++++++++++++++
        '''
        # estimate EW_fit_err
        # -------------------
        if num_refit > 0:
            # {'mod_gauss_amplitude',
            #  'mod_gauss_center',
            #  'mod_gauss_fwhm',
            #  'mod_gauss_sigma'}
            out_gauss_refit_amplitude = np.zeros(num_refit)
            out_gauss_refit_center = np.zeros(num_refit)
            out_gauss_refit_sigma = np.zeros(num_refit)
            # noise_fit = np.random.randn(num_refit,resi_range.size)*noise_std
            for i in xrange(int(num_refit)):
                # resi_range_with_noise = resi_range + noise_fit[i,:]
                resi_range_with_noise = resi_range + \
                                        np.random.randn(resi_range.size) * noise_std
                out_gauss_refit = mod_gauss.fit(resi_range_with_noise,
                                                par_gauss,
                                                x=wave_range)
                out_gauss_refit_amplitude[i],\
                out_gauss_refit_center[i],\
                out_gauss_refit_sigma[i] =\
                    out_gauss_refit.params[mod_gauss.prefix + 'amplitude'].value,\
                    out_gauss_refit.params[mod_gauss.prefix + 'center'].value,\
                    out_gauss_refit.params[mod_gauss.prefix + 'sigma'].value
            line_indx.update({'mod_gauss_amplitude_std': np.std(out_gauss_refit_amplitude),
                              'mod_gauss_center_std':    np.std(out_gauss_refit_center),
                              'mod_gauss_sigma_std':     np.std(out_gauss_refit_sigma)})
        return line_indx
    except Exception:
        return measure_line_index_null_result()  # necessary


def measure_line_index_null_result():
    """generate default value (nan/False) when measurement fails

    :return: default value (nan/False)
    """
    return {'SN_local_flux_err':        np.nan,
            'SN_local_flux_std':        np.nan,
            'flag_good_pixel':          False,
            'EW_int':                   np.nan,
            'EW_int_err':               np.nan,
            'mod_linear_slope':         np.nan,
            'mod_linear_slope_err':     np.nan,
            'mod_linear_intercept':     np.nan,
            'mod_linear_intercept_err': np.nan,
            'mod_gauss_amplitude':      np.nan,
            'mod_gauss_amplitude_err':  np.nan,
            'mod_gauss_center':         np.nan,
            'mod_gauss_center_err':     np.nan,
            'mod_gauss_sigma':          np.nan,
            'mod_gauss_sigma_err':      np.nan,
            'mod_gauss_amplitude_std':  np.nan,
            'mod_gauss_center_std':     np.nan,
            'mod_gauss_sigma_std':      np.nan}


# pure fit:    100 loops, best of 3: 8.06 ms per loop (1 int + 1 fit)
# 1000 re-fit: 1 loops, best of 3: 378 ms per loop (1 int + 1 fit + 100 re-fit)


def measure_line_index_loopfun(filepath):
    """loopfun for measuring line index

    :param filepath:    input file path
    :return:
    """
    num_refit = 50
    line_info_dib5780 = {'line_center':         5780,
                         'line_range':          (5775, 5785),
                         'line_shoulder_left':  (5755, 5775),
                         'line_shoulder_right': (5805, 5825)}
    line_info_dib5797 = {'line_center':         5797,
                         'line_range':          (5792, 5802),
                         'line_shoulder_left':  (5755, 5775),
                         'line_shoulder_right': (5805, 5825)}
    line_info_dib6284 = {'line_center':         6285,
                         'line_range':          (6280, 6290),
                         'line_shoulder_left':  (6260, 6280),
                         'line_shoulder_right': (6310, 6330)}
    try:
        # read spectrum
        # -------------
        spec = read_spectrum(filepath, 'auto')

        # measure DIBs
        # ------------
        # DIB5780
        line_indx_dib5780 = measure_line_index(wave=spec['wave'],
                                               flux=spec['flux'],
                                               flux_err=spec['flux_err'],
                                               mask=spec['and_mask'],
                                               line_info=line_info_dib5780,
                                               num_refit=num_refit,
                                               z=0)
        # DIB5797
        line_indx_dib5797 = measure_line_index(wave=spec['wave'],
                                               flux=spec['flux'],
                                               flux_err=spec['flux_err'],
                                               mask=spec['and_mask'],
                                               line_info=line_info_dib5797,
                                               num_refit=num_refit,
                                               z=0)
        # DIB6284 
        line_indx_dib6284 = measure_line_index(wave=spec['wave'],
                                               flux=spec['flux'],
                                               flux_err=spec['flux_err'],
                                               mask=spec['and_mask'],
                                               line_info=line_info_dib6284,
                                               num_refit=num_refit,
                                               z=0)
        return line_indx_dib5780, line_indx_dib5797, line_indx_dib6284
    except Exception:
        return measure_line_index_null_result(),\
               measure_line_index_null_result(),\
               measure_line_index_null_result()


def measure_line_index_recover_spectrum(wave, params, norm=False):
    from lmfit.models import LinearModel, GaussianModel
    mod_linear = LinearModel(prefix='mod_linear_')
    mod_gauss = GaussianModel(prefix='mod_gauss_')
    par_linear = mod_linear.make_params()
    par_gauss = mod_gauss.make_params()
    par_linear['mod_linear_slope'].value = params[0]
    par_linear['mod_linear_intercept'].value = params[1]
    par_gauss['mod_gauss_amplitude'].value = params[2]
    par_gauss['mod_gauss_center'].value = params[3]
    par_gauss['mod_gauss_sigma'].value = params[4]
    if not norm:
        flux = 1 - mod_gauss.eval(params=par_gauss, x=wave)
    else:
        flux = \
            (1 - mod_gauss.eval(params=par_gauss, x=wave)) * \
            mod_linear.eval(params=par_linear, x=wave)
    return flux


# %% test
if __name__ == '__main__':
    # filepath = 'spec-6064-56097-0980.fits'
    # filepath = 'spec-1230-52672-0233.fits'
    # filepath = 'spec-56309-GAC088N20V1_sp08-126.fits'
    # spec = read_spectrum(filepath, filesource)

    # fig = plt.figure()
    # ax = fig.add_axes()
    # plt.plot(spec['wave']/(1+z), spec['flux'], 'b')
    # plt.plot(spec['wave']/(1+z), spec['flux_err'], 'r')

    # from spectrum import generate_lamost_filepath,read_spectrum,measure_line_index
    filepath = '/home/cham/data/ToolFun/spectra/spec-1230-52672-0233.fits'
    filesource = 'auto'
    #    filepath = r'/pool/LAMOST/DR2/spectra/fits/F5902/spec-55859-F5902_sp01-001.fits'
    #    filesource = 'lamost_dr2'
    spec = read_spectrum(filepath, filesource)  # 10 loops, best of 3: 35.7 ms per loop
    #    line_indx_pack = measure_line_index_loopfun(filepath)
    z = 0.00205785
    line_info_dib6284 = {'line_center':         6285,
                         'line_range':          (6280, 6290),
                         'line_shoulder_left':  (6260, 6280),
                         'line_shoulder_right': (6310, 6330)}
    line_indx = measure_line_index(wave=spec['wave'],
                                   flux=spec['flux'],
                                   flux_err=spec['flux_err'],
                                   mask=spec['and_mask'],
                                   line_info=line_info_dib6284,
                                   num_refit=100,
                                   z=z)
    print line_indx
    '''
    45 ms for integration and other procedures
    380 ms for 100 refits
    In the fastest way (45ms), run 40 line indices on 4 million spectra:
    0.045*40*4E6/24/86400 ~ 3.5 days
    In the slowest way (380ms)
    0.420*40*4E6/24/86400 ~ 32.5 days
    '''

# %%
# plt.plot(wave_range, flux_range, 'bo-')
# plt.plot(wave_range, (1-mod_gauss.eval(par_gauss, x=wave_range))*cont_range,'k--')
# plt.plot(wave_range, (1-out_gauss.best_fit)*cont_range, 'r-')
# plt.show()

# %%
# plt.plot(wave_shoulder, flux_shoulder, 'bo-')
# plt.plot(wave_shoulder, mod_linear.eval(par_linear, x=wave_shoulder),'k--')
# plt.plot(wave_shoulder, out_linear.best_fit, 'r-')
# plt.show()

# %%
# fig = plt.figure()
# ax = fig.add_axes()
# plt.plot(\
#        wave,flux,'b-')
# plt.plot(\
#        wave_shoulder_left,\
#        flux_shoulder_left,\
#        'g-')
# plt.plot(\
#        wave_shoulder_right,\
#        flux_shoulder_right,\
#        'g-')
# plt.plot(\
#        hstack([wave_shoulder_left,wave_shoulder_right]),\
#        pget.best_fit,\
#        'r-')
# plt.plot(\
#        wave_range,flux_range,'b-')
#
# plt.xlim(6260,6330)
# plt.ylim(50,80)




# %%
# linear_fun_slope = \
#        (flux_shoulder_right_med-flux_shoulder_left_med)/\
#        (wave_shoulder_right_med-wave_shoulder_left_med)
# linear_fun_ofst  = \
#        flux_shoulder_left_med-linear_fun_slope*wave_shoulder_left_med
# linear_fun = lambda x: linear_fun_slope*x+linear_fun_ofst

# fig = plt.figure()
# ax = fig.add_axes()
# plt.plot(spec['wave']/(1+z), spec['flux'], 'b')
# plt.plot(spec['wave']/(1+z), spec['flux_err'], 'r')
# plt.plot(wave, linear_fun(wave), 'r')
# plt.plot(wave_shoulder_left, flux_shoulder_left, 'go')
# plt.plot(wave_shoulder_right, flux_shoulder_right, 'go')

# %% document

# % data structure +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# filesource:   'lamost_dr2' ----------------------------------------------
# http://dr2.lamost.org/doc/data-production-description#toc_3
#
#  RowNumber 	Data                Type
#  1           Flux                float
#  2           Inverse Variance 	float
#  3           WaveLength          float
#  4           Andmask             float
#  5           Ormask              float
#
#
# filesource:   'sdss_dr10' -----------------------------------------------
# http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
#
#  HDU 0  : Header info from spPlate
#  HDU 1  : Coadded spectrum from spPlate --> use this
#  HDU 2  : Summary metadata copied from spAll
#  HDU 3  : Line fitting metadata from spZline
#  HDU 4+ : [Optional] Individual spCFrame spectra [B, R for each exposure]
#
#
# HDU 0: Header keywords only
#
# Copied from spPlate with the following additions/modifications:
#
#    PLUG_RA, PLUG_DEC, THING_ID, FIBERID: added from spAll
#    NEXP and EXPID*: modified to just include the frames which contributed to this fiber
#    Removed keywords which apply only to single exposures
#
# HDU 1 (extname COADD): Coadded Spectrum from spPlate
#
# Binary table with columns:
#  Required    Columns
#  Col     Name        Type        Comment
#  1       flux        float32 	coadded calibrated flux [10-17 ergs/s/cm2/Å]
#  2       loglam      float32 	log10(wavelength [Å])
#  3       ivar        float32 	inverse variance of flux
#  4       and_mask 	int32       AND mask
#  5       or_mask 	int32       OR mask
#  6       wdisp       float32 	wavelength dispersion in pixel=dloglam units
#  7       sky         float32 	subtracted sky flux [10-17 ergs/s/cm2/Å]
#  8       model       float32 	pipeline best model fit used for classification and redshift
