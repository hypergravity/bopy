# -*- coding: utf-8 -*-
"""

Author
------
Bo Zhang, migrated module

Email
-----
bozhang@nao.cas.cn

Created on
----------
- Fri Jun 15 13:16:00 2016    deredden

Modifications
-------------
- Fri Jun 15 13:16:00 2016    re-format code

Aims
----
- to de-redden the extra-galactic spectra using CCM 1989/1994 law
  http://adsabs.harvard.edu/abs/1989ApJ...345..245C
  lambda range; 0.125 um (1250 A) < lambda < 3.5 um (35000 A)

"""

import numpy as np


def ccm(wave, Av=1., Rv=3.1, strict_ccm=False):
    """ Returns the Cardelli, Clayton, and Mathis (CCM) reddening curve.

    Parameters
    ----------
    wave: numpy.ndarray, float
        the wavelength array, in \AA unit

    Returns
    -------

    """
    # Convert to inverse microns (um^-1)
    x = 10000. / wave
    a = 0.0 * x
    b = 0.0 * x

    # 1. Infrared
    good = np.greater_equal(x, 0.3) * np.less(x, 1.1)
    if len(np.nonzero(good)) > 0:
        a = np.where(good,  0.574 * np.power(x, 1.61), a)
        b = np.where(good, -0.527 * np.power(x, 1.61), b)

    # 2. Optical/NIR
    good = np.greater_equal(x, 1.1) * np.less(x, 3.3)
    if len(np.nonzero(good)) > 0:  # Use new constants from O'Donnell (1994)
        y = x - 1.82
        if strict_ccm:
            # Original coefficients from CCM 1989
            c1 = [0.32999, -0.7753, 0.01979, 0.72085, -0.02427,
                  -0.50447, 0.17699, 1.0]
            c2 = [-2.09002, 5.3026, -0.62251, -5.38434, 1.07233,
                  2.28305, 1.41338, 0.0]
        else:
            # New coefficents from O'Donnell (1994)
            c1 = [-0.505, 1.647, -0.827, -1.718, 1.137, 0.701,
                  -0.609, 0.104, 1.0]
            c2 = [3.347, -10.805, 5.491, 11.102, -7.985, -3.989,
                  2.908, 1.952, 0.0]

        a = np.where(good, np.polyval(c1, y), a)
        b = np.where(good, np.polyval(c2, y), b)

    # 3. Mid-UV
    good = np.greater_equal(x, 3.3) * np.less(x, 8.)
    if len(np.nonzero(good)) > 0:
        good1 = np.greater(x, 5.9)

        F_a = np.polyval([-0.009779, -0.04473, 0., 0.], x - 5.9)
        F_b = np.polyval([0.1207, 0.2130, 0., 0.], x - 5.9)
        F_a = np.where(good1, F_a, 0.0)
        F_b = np.where(good1, F_b, 0.0)

        a_ = 1.752 - 0.316 * x - (0.104 / (np.power(x - 4.67, 2.) + 0.341))
        b_ = -3.090 + 1.825 * x + (1.206 / (np.power(x - 4.67, 2) + 0.263))
        a = np.where(good, a_ + F_a, a)
        b = np.where(good, b_ + F_b, b)

    # 4. Far-UV
    good = np.greater_equal(x, 8) * np.less_equal(x, 11)
    if len(np.nonzero(good)) > 0:
        y = x - 8.0
        c1 = [-0.07, 0.137, -0.628, -1.073]
        c2 = [0.374, -0.42, 4.257, 13.67]
        a = np.where(good, np.polyval(c1, y), a)
        b = np.where(good, np.polyval(c2, y), b)

    # A(lambda)/A(V) = a(x) + b(x)/Rv, see CCM 1989 eq. (1)
    Al = Av * (a + b / Rv)
    return Al


def dered_ccm(wave, flux, flux_err=None, Av=1.,
              Rv=3.1, z=None, strict_ccm=True):
    """ un-redden the spectra by specified fore-groud extinction

    Parameters
    ----------
    wave: array_like
        wavelength array

    flux: array_like
        flux array

    Av: float
        the fore-ground extinction

    Rv: float
        the Rv parameters in CCM extinction law

    strict_ccm: bool
        If True, use CCM(1989), else use O'Donnell(1994) for optical band.

    Returns
    -------
    flux_unred: array_like
        the un-reddened flux array


    --------------------------------------------------------------
    This function is migrated from IDL routine,
    here below goes the original doc-sting:
    --------------------------------------------------------------

    pro ccm_UNRED, wave, flux, ebv, funred, R_V = r_v
   
    NAME:
        CCM_UNRED
    PURPOSE:
        Deredden a flux vector using the CCM 1989 parameterization
    EXPLANATION:
        The reddening curve is that of Cardelli, Clayton, and Mathis (1989 ApJ.
        345, 245), including the update for the near-UV given by O'Donnell
        (1994, ApJ, 422, 158).   Parameterization is valid from the IR to the
        far-UV (3.5 microns to 0.1 microns).

        Users might wish to consider using the alternate procedure FM_UNRED
        which uses the extinction curve of Fitzpatrick (1999).
    CALLING SEQUENCE:
        CCM_UNRED, wave, flux, ebv, funred, [ R_V = ]
               or
        CCM_UNRED, wave, flux, ebv, [ R_V = ]
    INPUT:
        WAVE - wavelength vector (Angstroms)
        FLUX - calibrated flux vector, same number of elements as WAVE
               If only 3 parameters are supplied, then this vector will
               updated on output to contain the dereddened flux.
        EBV  - color excess E(B-V), scalar.  If a negative EBV is supplied,
               then fluxes will be reddened rather than deredenned.

    OUTPUT:
        FUNRED - unreddened flux vector, same units and number of elements
               as FLUX

    OPTIONAL INPUT KEYWORD
        R_V - scalar specifying the ratio of total selective extinction
               R(V) = A(V) / E(B - V).    If not specified, then R_V = 3.1
               Extreme values of R(V) range from 2.75 to 5.3

    EXAMPLE:
        Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
        is altered by a reddening of E(B-V) = 0.1.   Assume an "average"
        reddening for the diffuse interstellar medium (R(V) = 3.1)

        IDL> w = 1200 + findgen(40)*50      ;Create a wavelength vector
        IDL> f = w*0 + 1                    ;Create a "flat" flux vector
        IDL> ccm_unred, w, f, -0.1, fnew  ;Redden (negative E(B-V)) flux vector
        IDL> plot,w,fnew

    NOTES:
        (1) The CCM curve shows good agreement with the Savage & Mathis (1979)
               ultraviolet curve shortward of 1400 A, but is probably
               preferable between 1200 and 1400 A.
        (2)  Many sightlines with peculiar ultraviolet interstellar extinction
               can be represented with a CCM curve, if the proper value of
               R(V) is supplied.
        (3)  Curve is extrapolated between 912 and 1000 A as suggested by
               Longo et al. (1989, ApJ, 339,474)
        (4) Use the 4 parameter calling sequence if you wish to save the
                 original flux vector.
        (5) Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
               curve (3.3 -- 8.0 um-1).    But since their revised curve does
               not connect smoothly with longer and shorter wavelengths, it is
               not included here.

    REVISION HISTORY:
        Written   W. Landsman        Hughes/STX   January, 1992
        Extrapolate curve for wavelengths between 900 and 1000 A   Dec. 1993
        Use updated coefficients for near-UV from O'Donnell   Feb 1994
        Allow 3 parameter calling sequence      April 1998
        Converted to IDLV5.0                    April 1998

    Added a redshift argument so that we can (de)redden in a different rest
    frame. It is used in the sense that the wavelengths are redshifted by z
    before the calculation is performed.
    We do this because the SED templates are defined in the observed frame.
    """


    if z is not None:
        wave = wave * (1 + z)

    # evaluate A_lambda from CCM extinction law
    Al = ccm(wave, Av=Av, Rv=Rv, strict_ccm=strict_ccm)

    # Now apply extinction correction to the input flux vector
    flux_dered = flux * np.power(10.0, 0.4 * Al)

    if flux_err is not None:
        # Now apply extinction correction to the input flux_err vector
        flux_err_dered = flux_err * np.power(10.0, 0.4 * Al)
        return flux_dered, flux_err_dered

    return flux_dered


# def R_z(wave, z, R_V=3.1, strict_ccm=0):
#     '''Compute host reddening law for -effective observed wavelength wave, at redshift
#    z, assuming R_V.'''
#     wave = wave / (1 + z)
#     a, b = ccm(wave, strict_ccm)
#     R = R_V * a + b
#     return (R, a, b)
