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
# todo : norm spectra --> plot spectra --> print file

def spec_quick_view(ax, spec_chunk_list, ranges=None,
                    num_spec_perpage=30, offset_perspec=0.5, wave_label='abs',
                    cl='', *args, **kwargs):

    return 0


def test_spec_quick_view():
    fig = plt.figure('test spec_quick_view method')
    # load a list of spectra

    # break into chunks

    # plot them

    print('@Cham: test spec_quick_view OK!')
    return 0


if __name__ == '__main__':
    test_spec_quick_view()

