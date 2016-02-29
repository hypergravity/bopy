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
from .continuum_normalization import (_cont_norm_running_quantile,
                                      _cont_norm_running_quantile_regions)


class Spec(Table):

    def __init__(self, *args, **kwargs):
        super(Spec, self).__init__(*args, **kwargs)
        assert 'wave' in self.colnames
        assert 'flux' in self.colnames

    def norm_spec(self, ranges=None, q=0.90, delta_lambda=100.,
                  norm_flux_colname='flux_norm'):
        flux_norm, ivar_norm = \
            continuum_normalize_training_q(self, q=0.90, delta_lambda=100.)
        self.add_columns([Column(flux_norm, 'flux_norm'),
                          Column(ivar_norm, 'ivar_norm')])


# ############################# #
# spectrum quick initialization #
# ############################# #

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


# ###################################### #
# continuum normalization for a spectrum #
# ###################################### #

def wave2ranges(wave, wave_intervals):
    """ convert wavelength intervals to (pixel) ranges """
    wave_intervals = np.array(wave_intervals)

    # assert wave_intervals is a 2-column array
    assert wave_intervals.shape[1] == 2

    ranges = np.zeros_like(wave_intervals)
    for i in xrange(len(wave_intervals)):
        ranges[i, 0] = np.sum(wave < wave_intervals[i, 0])
        ranges[i, 1] = np.sum(wave < wave_intervals[i, 2]) - 1
    return ranges


def continuum_normalize_training_q(spec, ranges=None, q=0.90, delta_lambda=100.):
        """ Continuum normalize the training set using a running quantile
            migrated from TheCannon

        Parameters
        ----------
        q: float
            The quantile cut (q between 0.0 and 1.0)
        delta_lambda: float
            The width of the pixel range used to calculate the median

        Returns
        -------
        flux_norm, ivar_norm

        """
        # print("Continuum normalizing the tr set using running quantile...")
        print('@Cham:: normalizing spectra using running quantile ...')
        if 'ivar' not in spec.colnames:
            # this is a bear spectrum without ivar data
            # produce an all-one array to replace the ivar
            spec.add_column(Column(np.ones_like(spec['flux']), 'ivar'))

        if ranges is None:
            # continuous spectrum
            return _cont_norm_running_quantile(
                spec['wave'], spec['flux'], spec['ivar'],
                q=q, delta_lambda=delta_lambda)
        else:
            return _cont_norm_running_quantile_regions(
                spec['wave'], spec['flux'], spec['ivar'],
                q=q, delta_lambda=delta_lambda, ranges=ranges)


if __name__ == '__main__':
    _test_spec_quick_init()