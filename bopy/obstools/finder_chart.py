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
- Sat Jan 16 19:55:57 2016

Modifications
-------------
- Sat Jan 16 19:55:57 2016    add finder_chart.py

Aims
----
- auto-generate stsci DSS finder chart

"""


import numpy as np
from astropy.coordinates import SkyCoord


def finder_chart_url_stscidss(
        ra=62.165336,
        dec=47.712509,
        v='poss2ukstu_red',
        epoch='J2000',
        h=5.,
        w=5.,
        f='gif'):
    """ return the url of stsci DSS image

    Parameters
    ----------

    ra: number/array
        right ascension. range: [0.0 360.0]

    dec: number/array.
        declination. range: [-90.0 +90.0]

    v: string
        {
        'poss2ukstu_red', 'poss2ukstu_blue', 'poss2ukstu_ir'
        'poss1_red', 'poss1_blue',
        'quickv',
        'phase2_gsc2', 'phase2_gsc1'
        }

    epoch: string {'J2000' 'B1950'}
        epoch of coordinate system

    h,w: number
        height and width of image. maximum is 60. (arcmin)

    f: string {'gif' 'fits'}
        file format

    Returns
    -------

    urlstr: string
        the link of the image
    """

    ra, dec = np.array(ra), np.array(dec)
    c = SkyCoord(ra, dec, unit="deg")

    urlform = ['https://archive.stsci.edu/cgi-bin/dss_search?v=%s',
               '&r=%02d+%02d+%f',
               '&d=%s%02d+%02d+%f',
               '&e=%s',
               '&h=%f',
               '&w=%f',
               '&f=%s',
               '&c=none&fov=NONE&v3=']
    urldict = {True: '%2B', False: '-'}

    if np.isscalar(c.ra.hms.h):
        # only 1 pair of [ra, dec]
        urlpiece = [
            urlform[0] % v,
            urlform[1] % (c.ra.hms.h, c.ra.hms.m, c.ra.hms.s),
            urlform[2] % (urldict[c.dec.dms.d >= 0],
                          np.abs(c.dec.dms.d),
                          np.abs(c.dec.dms.m),
                          np.abs(c.dec.dms.s)),
            urlform[3] % epoch,
            urlform[4] % h,
            urlform[5] % w,
            urlform[6] % f,
            urlform[7]
        ]
        urlstr = ''.join(urlpiece)
        print c.dec.dms.m
        return urlstr
    else:
        # a list of pairs of [ra, dec]
        urlstr = list()
        for c_ in c:
            urlpiece = [
                urlform[0] % v,
                urlform[1] % (c_.ra.hms.h, c_.ra.hms.m, c_.ra.hms.s),
                urlform[2] % (urldict[c_.dec.dms.d >= 0],
                              np.abs(c_.dec.dms.d),
                              np.abs(c_.dec.dms.m),
                              np.abs(c_.dec.dms.s)),
                urlform[3] % epoch,
                urlform[4] % h,
                urlform[5] % w,
                urlform[6] % f,
                urlform[7]
            ]
            urlstr.append(''.join(urlpiece))
        urlstr = np.array(urlstr)
    return urlstr


if __name__ == '__main__':
    # test finder_chart
    print finder_chart_url_stscidss()
    print '@Cham: test passed ...'
