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
- Tue Mar  8 15:26:00 2016    read_spectrum

Modifications
-------------
-

Aims
----
- read various kinds of spectra

"""

import os
import numpy as np
from astropy.io import fits
from .spec import Spec


def reconstruct_wcs_coord_from_fits_header(hdr, dim=1):
    """ reconstruct wcs coordinates (e.g., wavelenght array) """
    # assert dim is not larger than limit
    assert dim <= hdr['NAXIS']

    # get keywords
    crval = hdr['CRVAL%d' % dim]
    cdelt = hdr['CDELT%d' % dim]
    crpix = hdr['CRPIX%d' % dim]
    naxis = hdr['NAXIS%d' % dim]

    # reconstruct wcs coordinates
    coord = np.arange(1 - crpix, naxis + 1 - crpix) * cdelt + crval
    return coord


def read_spectrum_elodie_r42000(fp):
    """ read spectrum from ELODIE library (R42000) """
    # assert the file exists
    assert os.path.exists(fp)

    # read fits
    hl = fits.open(fp)

    # reconstruct wave array
    wave = reconstruct_wcs_coord_from_fits_header(hl[0].header, dim=1)
    # flux
    flux = hl[0].data
    # flux err
    flux_err = hl[2].data
    # flux ivar
    flux_ivar = 1 / flux_err ** 2.

    # reconstruct spec
    sp = Spec(data=[wave, flux, flux_ivar, flux_err],
              names=['wave', 'flux', 'flux_ivar', 'flux_err'])
    return sp
