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
- Sat Jan 23 11:38:45 2016

Modifications
-------------
-

Aims
----
- using pyraf.iraf.onedspec.wspectext routine to write spectra in text format

"""

import os
import sys
from pyraf import iraf
from astropy.io import fits


def wspectext_allap(filepath, apaxis=0):
    """ this can be directly evaluated in terminal

    Parameters
    ----------
    filepath: string
        filepath

    apaxis: integer
        0 or 1.
        0 denotes dispersion along rows (this is the case for Xinglong216HRS)
        1 denotes dispersion along columns

    Examples
    --------
    wspectext_allap('./w20160120001s.fits')

    """
    nrow, ncol = fits.open(filepath)[0].data.shape
    print '@Cham: nrow = %d, ncol=%d' % (nrow, ncol)
    
    iraf.noao()
    iraf.noao.onedspec()
    
    filename = filepath.split(os.path.sep)[-1]
    
    # determine dirname & dirpath
    dirname = filepath.split(os.path.sep)[-1].replace('.fits', '')
    if os.path.dirname(filepath) == '':
        dirpath = dirname
    else:
        dirpath = os.path.dirname(filepath) + os.path.sep + dirname

    if os.path.exists(dirpath):
        # if dirpath exists
        print '@Cham: directory exists ... (%s)' % dirpath
    else:
        # if dirpath doesn't exist, mkdir
        os.mkdir(dirpath)
        print '@Cham: mkdir %s'%dirpath
    
    # execute wspectest
    if not apaxis == 1:
        for apnum in xrange(1, nrow+1):
            _input = '%s[1:%d, %d]' % (filename, ncol, apnum)
            _output = dirpath + os.sep + filename.replace('.fits', '_%04d.dat' % apnum)
            print '@Cham: wspectext running ... [%s]' % _output
            iraf.onedspec.wspectext(input=_input, output=_output, header='no')
    else:
        for apnum in xrange(1, ncol+1):
            _input = '%s[%d, 1:%d]' % (filename, apnum, nrow)
            _output = dirpath + os.sep + filename.replace('.fits', '_%04d.dat' % apnum)
            print '@Cham: wspectext running ... [%s]' % _output
            iraf.onedspec.wspectext(input=_input, output=_output, header='no')

    # print masseges
    print '@Cham: ----------------------'
    print '@Cham: mission complete!'
    print '@Cham: ----------------------'


if __name__ == '__main__':
    """ this can be directly evaluated in terminal

    Examples
    --------
    ipython py_wspectxt ./w*s.fits
    or:
    python py_wspectxt ./w*s.fits 0

    """

    if len(sys.argv) < 2:
        print '@Cham: not enough arguments ...'
    elif len(sys.argv) == 2:
        wspectext_allap(sys.argv[1])
    elif len(sys.argv) >= 3:
        if sys.argv[-1]=='0' or sys.argv[-1]=='1':
            # apaxis specified
            for i in xrange(1, len(sys.argv)-1):
                wspectext_allap(sys.argv[i], sys.argv[-1])
        else:
            # apaxis not specified
            for i in xrange(1, len(sys.argv)):
                wspectext_allap(sys.argv[i])
