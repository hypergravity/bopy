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
- Sun Feb  28 15:02:06 2016    Spec

Modifications
-------------
- Fri Nov 20 10:16:59 2015    reformatting code

Aims
----
- spec class

"""


import numpy as np
from astropy.io import fits
from astropy.table import Table, Column


class Spec(Table):

    def __init__(self, *args, **kwargs):
        super(Spec, self).__init__(*args, **kwargs)
        assert 'wave' in self.colnames
        assert 'flux' in self.colnames

    # TODO: implement several functions, e.g.,
    def norm_flux(self, q, norm_flux_colname='flux_norm'):
        # self[norm_flux_colname] =
        return 0


def spec_quick_init(wave, flux):
    """

    Parameters
    ----------
    wave : numpy.ndarray
        wavelength array

    flux : numpy.ndarray
        flux array

    Returns
    -------
    spec_ : bopy.spec.Spec
        Spec, at least contains 'wave' and 'flux' columns.

    """
    spec_ = Spec(
        [Column(np.array(wave), 'wave'),
         Column(np.array(flux), 'flux')])
    return spec_


def _test_spec_quick_init():
    wave = np.arange(5000., 7000., 0.91)
    flux = np.sin(wave/1000)*1000 + 30000. + np.random.rand(len(wave))
    spec_qi = spec_quick_init(wave, flux)
    spec_qi.pprint()
    print '--------------------------------------'
    print '@Cham: _test_spec_quick_init() OK ...'
    print '--------------------------------------'

if __name__ == '__main__':
    _test_spec_quick_init()