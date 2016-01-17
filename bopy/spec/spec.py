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


import numpy as np
from astropy.io import fits
from astropy.table import Table, Column



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


def _test_lamost_filepath():
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

    # LAMOST DR2/DR3 spectrum
    if filesource == 'lamost_dr3':
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


def _test_read_spectrum():
    fp = '/home/cham/PycharmProjects/bopy/bopy/data/test_spectra/lamost_dr3/'\
         + lamost_filepath('GAC_061N46_V3', 55939, 7, 78)
    print fp
    sp = read_spectrum(fp)
    sp.pprint()


# -------------------------------------------------------------------------
# Test the module ...
if __name__ == '__main__':
    print ''
    print '@Cham: start to test the module ...'
    print ''
    print '@Cham: testing ''lamost_filepath'' ...'
    _test_lamost_filepath()
    print ''
    print '@Cham: testing ''read_spectrum'' ...'
    _test_read_spectrum()
    print '@Cham: OK'


# documents of data structures (LAMOST and SDSS spectra)
#
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
