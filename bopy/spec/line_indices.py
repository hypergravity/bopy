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
- Fri Jul  3 13:13:06 2015    read_spectrum

Modifications
-------------
- Wed Jul 29 21:46:00 2015    measure_line_index
- Fri Nov 20 10:16:59 2015    reformatting code
- Sat Jan 16 19:55:57 2016    migrate from spec.py
- Thu Jul 14 23:57:57 2016    plot every line indice

Aims
----
- measure line index from spectra

"""

import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import LinearModel, GaussianModel
from .lamost import read_spectrum

def measure_line_index(filepath, wave, flux, flux_err=None, mask=None,
                       line_info=None, num_refit=0, z=0.):
    """Measure line index / line EW and have it plotted

    Parameters
    ----------

    filepath: string
        path of the spec document
        
    wave: array
        wavelength vector

    flux: array
        flux vector

    flux_err: array
        flux error vector (optional)
        If un-specified, auto-generate an np.ones array

    mask: array
        andmask or ormask (optional)
        If un-specified, auto-generate an np.ones array (evenly weighted)

    line_info: dict
        information about spectral line, eg:
        line_info_dib5780 = {'line_center':         5780,
                             'line_range':          (5775, 5785),
                             'line_shoulder_left':  (5755, 5775),
                             'line_shoulder_right': (5805, 5825)}

    num_refit: non-negative integer
        number of refitting.
        If 0, no refit will be performed
        If positive, refits will be performed after adding random noise

    z: float
        redshift (only specify when z is large)

    Returns
    -------

    line_indx: dict
        A dictionary type result of line index.
        If any problem encountered, return the default result (filled with nan).

    """ 
    try:
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
        if not z == 0:
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
        ''' 5. estimate continuum ---------------------------------------------
        '''
        cont_shoulder = out_linear.best_fit
        noise_std = np.std(flux_shoulder / cont_shoulder)
        cont_range = mod_linear.eval(out_linear.params, x=wave_range)
        resi_range = 1 - flux_range / cont_range

        ''' 6.1 Integrated EW -------------------------------------------------
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

        ''' 6.2 gaussian model ------------------------------------------------
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

        ''' 7. add noise and re-fit -------------------------------------------
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

        ''' 8. plot and save image --------------------------------------------
        '''
        save_image_line_indice(filepath, wave, flux, ind_range, cont_range, 
                               ind_shoulder, line_info)
        return line_indx
    except Exception:
        return measure_line_index_null_result()  # necessary


def save_image_line_indice(filepath, wave, flux, ind_range, cont_range, 
                           ind_shoulder, line_info):
    """Plot a line indice and save it as a .png document.
    
    Parameters
    ----------
    filepath: string
        path of the spec document
    
    wave: array
        wavelength vector
        
    flux: array
        flux vector
        
    ind_range: array
        bool indicating the middle range of a particular line
        
    cont_range: array
        continuum flux of the middle range derived from linear model
        
    ind_shoulder: array
        bool indicating the shoulder range of a particular line
    
    line_info: dict
        information about spectral line, eg:
        line_info_dib5780 = {'line_center':         5780,
                             'line_range':          (5775, 5785),
                             'line_shoulder_left':  (5755, 5775),
                             'line_shoulder_right': (5805, 5825)}

    """
    #   We suppose that the filepath is 'Users/bo/spec-11111-GAC_082N27_B1_sp01-001.fits'
    #   Then filepath[9:-5] = 'spec-11111-GAC_082N27_B1_sp01-001'
    #   We want to save it in 'Users/bo/images'
    #   name of the document is 'spec-11111-GAC_082N27_B1_sp01-001_line5780.png'
    fig = plt.figure()
    plt.plot(wave[ind_range], flux[ind_range], 'r-')
    plt.plot(wave[ind_range], cont_range, 'b-')
    plt.plot(wave[ind_shoulder],flux[ind_shoulder], 'm-')
    plt.title(r'line' + str(line_info['line_center']) + r'of ' + filepath[9:-5])
    add = filepath[9:-5] + '_line' + str(line_info['line_center'])
    plt.title(r'line' + str(line_info['line_center']) + r'of ' + filepath[31:-1])
    # fig.show()
    fig.savefig('/Users/bo/images/'+add+'.png')


def measure_line_index_null_result():
    """generate default value (nan/False) when measurement fails

    Returns
    -------
    default value (nan/False)
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
    Parameters
    ----------
    filepath: string
        path of the spec document
    
    Returns
    -------
    several line_indx: tuple
        every line_indx is a dictionary type result of line index.
    """
    num_refit = 0
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


def get_spectra_path():
    """Get the path of several spectra.
    
    Returns
    -------
    filename: list
        filepathes of all the spectra in finder '/Users/bo/spec/'
    """
    rootdir = '/Users/bo/spec'
    filename = []
    for parent,dirnames,filenames in os.walk(rootdir):
        for filename in filenames:
            filename.append(os.path.join(parent,filename))
    n_plus1 = len(filename)
    filename = filename[1:n_plus1]
    return filename
    

def test_():
    filepath = get_spectra_path()
    filesource = 'auto'
    #    filepath = r'/pool/LAMOST/DR2/spectra/fits/F5902/spec-55859-F5902_sp01-001.fits'
    #    filesource = 'lamost_dr2'
    spec = read_spectrum(filepath[0], filesource)  # 10 loops, best of 3: 35.7 ms per loop
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


def test_measure_line_index():
    filepath = get_spectra_path()
    n = len(filepath)
    line_indx_star = [[]for i in range(3)]
    for i in range(n):
        line_indx = measure_line_index_loopfun(filepath[i])
        line_indx_star[0].append(line_indx[0])
        line_indx_star[1].append(line_indx[1])
        line_indx_star[2].append(line_indx[2])
    return line_indx_star


def get_equivalent_width(line_indx_star):
    EW = [[]for i in range(3)]
    n = len(line_indx_star[0])
    for i in range(3):
        for j in range(n):
            EW[i].append(line_indx_star[i][j]['EW_int'])
    return EW
    

def plot_equivalent_width_hist(EW_star):
    titles = ["5780","5797","6285"]
    fig, axes = plt.subplots(1, 3, figsize=(8, 8))
    for i in range(3):
        ax = axes[0, i]
        ax.hist(EW_star[i],facecolor='red',alpha=0.5)
        ax.set_xlabel('equivalent width')
        ax.set_ylabel('number')
        ax.set_title('Histogram of equivalent width_line'+titles[i])
        plt.tight_layout()
    plt.show()

    
def plot_line_indices(EW_star):
    titles = ["5780","5797","6285"]
    fig, axes = plt.subplots(3,3,figsize=(64, 64))
    for i in range(3):
        for j in range(i+1):
            ax = axes[i,j]
            ax.set_title(titles[i]+" - "+titles[j],fontsize = 8)
            ax.set_ylabel(titles[i], fontsize=8)
            ax.set_xlabel(titles[j], fontsize=8)
            ax.plot(EW_star[j],EW_star[i],'ob',markersize=3, alpha=0.5)
    plt.tight_layout()


# %% test
if __name__ == '__main__':
    line_indx_star = test_measure_line_index()
    EW_star = get_equivalent_width(line_indx_star)
    #plot_equivalent_width_hist(EW_star)
    #plot_line_indices(EW_star)
