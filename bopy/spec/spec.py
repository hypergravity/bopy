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
- Fri Nov 20 10:16:59 2015    reformatting code

Aims
----
- generate LAMOST spectra file name/path
- to read LAMOST/SDSS spectra

"""

from astropy.io import fits
from astropy.table import Table, Column
import numpy as np


def lamost_filepath(planid, mjd, spid, fiberid, dirpath=''):
    """ generate file path of A LAMOST spectrum

    Parameters
    ----------

    planid: string
        planid

    mjd: 5-digit integer
        mjd

    spid: 2-digit integer
        spid, the number of the spectrogragh

    fiberid: 4-digit integer
        fiberid

    dirpath: string
        the root directory for storing spectra.

    Rreturns
    --------

    filepath: string
        the path of root dir of directory (prefix).
        if un-specified, return file name.

    """

    if dirpath == '' or dirpath is None:
        # return file name
        if np.isscalar(mjd):
            # if only input one item
            return 'spec-%05d-%s_sp%02d-%03d.fits' % (mjd, planid, spid, fiberid)
        else:
            # if input a list of items
            return np.array(['spec-%05d-%s_sp%02d-%03d.fits' %
                             (mjd[i], planid[i], spid[i], fiberid[i])
                             for i in xrange(len(mjd))])
    else:
        # return file path
        if not dirpath[-1] == '/':
            dirpath += '/'

        if np.isscalar(mjd):
            # if only input one item
            return '%s%s/spec-%05d-%s_sp%02d-%03d.fits' % (dirpath, planid, mjd, planid, spid, fiberid)
        else:
            # if input a list of items
            return np.array(['%s%s/spec-%05d-%s_sp%02d-%03d.fits' %
                             (dirpath, planid[i], mjd[i], planid[i], spid[i], fiberid[i])
                             for i in xrange(len(mjd))])


def test_lamost_filepath():
    """test function **lamost_filepath**
    """
    print lamost_filepath('GAC_061N46_V3', 55939, 7, 78)
    print lamost_filepath('GAC_061N46_V3', 55939, 7, 78, '/')


def read_spectrum(filepath, filesource='auto'):
    """read SDSS/LAMOST spectrum

    Parameters
    ----------

    filepath: string
        input file path

    filesource: string
        {'sdss_dr12' / 'lamost_dr2' / 'lamost_dr3'}

    Returns
    -------

    specdata: astropy.table.Table
        spectra as a table

    """
    # auto-identify the spectrum origination
    if filesource == 'auto':
        telescope = fits.open(filepath)[0].header['TELESCOP']
        if telescope == 'SDSS 2.5-M':
            return read_spectrum(filepath, filesource='sdss_dr12')
        if telescope == 'LAMOST':
            return read_spectrum(filepath, filesource='lamost_dr3')
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
        # for flux_err, convert inf to nan
        flux_err[np.isinf(flux_err.data)] = np.nan
        return Table([wave, flux, flux_err, ivar, and_mask, or_mask])
    return None


def test_read_spectrum():
    fp = '/home/cham/PycharmProjects/bopy/bopy/data/test_spectra/lamost_dr3/'\
         + lamost_filepath('GAC_061N46_V3', 55939, 7, 78)
    print fp
    sp = read_spectrum(fp, filesource='lamost_dr2')
    sp.pprint()


# -------------------------------------------------------------------------
# Test the module ...
if __name__ == '__main__':
    print ''
    print '@Cham: start to test the module ...'
    print ''
    print '@Cham: testing ''lamost_filepath'' ...'
    test_lamost_filepath()
    print ''
    print '@Cham: testing ''read_spectrum'' ...'
    test_read_spectrum()

    # # filepath = 'spec-6064-56097-0980.fits'
    # # filepath = 'spec-1230-52672-0233.fits'
    # # filepath = 'spec-56309-GAC088N20V1_sp08-126.fits'
    # # spec = read_spectrum(filepath, filesource)
    #
    # # fig = plt.figure()
    # # ax = fig.add_axes()
    # # plt.plot(spec['wave']/(1+z), spec['flux'], 'b')
    # # plt.plot(spec['wave']/(1+z), spec['flux_err'], 'r')
    #
    # # from spectrum import generate_lamost_filepath,read_spectrum,measure_line_index
    # filepath = '/home/cham/data/ToolFun/spectra/spec-1230-52672-0233.fits'
    # filesource = 'auto'
    # #    filepath = r'/pool/LAMOST/DR2/spectra/fits/F5902/spec-55859-F5902_sp01-001.fits'
    # #    filesource = 'lamost_dr2'
    # spec = read_spectrum(filepath, filesource)  # 10 loops, best of 3: 35.7 ms per loop
    # #    line_indx_pack = measure_line_index_loopfun(filepath)
    # z = 0.00205785
    # line_info_dib6284 = {'line_center':         6285,
    #                      'line_range':          (6280, 6290),
    #                      'line_shoulder_left':  (6260, 6280),
    #                      'line_shoulder_right': (6310, 6330)}
    # line_indx = measure_line_index(wave=spec['wave'],
    #                                flux=spec['flux'],
    #                                flux_err=spec['flux_err'],
    #                                mask=spec['and_mask'],
    #                                line_info=line_info_dib6284,
    #                                num_refit=100,
    #                                z=z)
    # print line_indx
    # '''
    # 45 ms for integration and other procedures
    # 380 ms for 100 refits
    # In the fastest way (45ms), run 40 line indices on 4 million spectra:
    # 0.045*40*4E6/24/86400 ~ 3.5 days
    # In the slowest way (380ms)
    # 0.420*40*4E6/24/86400 ~ 32.5 days
    # '''

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
