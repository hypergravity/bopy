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
- Sun Feb 28 14:39:16 2016    migrated to bopy.spec.lamost
- Fri Jul 15 16:08:00 2016    migrate read_spectrum to read_spectrum.py


Aims
----
- generate LAMOST spectra file name/path

"""

# from __future__ import print_function
import numpy as np
# from astropy.io import fits
# from astropy.table import Table, Column


def lamost_filepath(planid, mjd, spid, fiberid, dirpath=''):
    """ generate file path of a LAMOST spectrum

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
    print (lamost_filepath('GAC_061N46_V3', 55939, 7, 78))
    print (lamost_filepath('GAC_061N46_V3', 55939, 7, 78, '/'))


# -------------------------------------------------------------------------
# Test the module ...
if __name__ == '__main__':
    print('')
    print('@Cham: start to test the module ...')
    print('')
    print('@Cham: testing ''lamost_filepath'' ...')
    _test_lamost_filepath()
    print('@Cham: OK')
