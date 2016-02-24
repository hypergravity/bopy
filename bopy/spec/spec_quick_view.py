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
- Wed Feb  24 15:05:00 2016     spec_quick_view

Modifications
-------------
-

Aims
----
- method spec_quick_view: tool for quick view of spectra

"""

import numpy as np
import matplotlib.pyplot as plt


# input: a list of spectra
# to do: norm spectra --> plot spectra --> print file

def spec_quick_view(spec_list, norm_flag='', *args):



def norm_spec(wave, flux, kind='median', amp=2.):
    """

    Args:
        wave:
        flux:
        kind: {median/med | max | }
        amp:

    Returns:

    """