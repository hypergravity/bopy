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
from astropy.table import Table
from .spec import spec_quick_init, Spec, norm_spec_pixel

# input: a list of spectra
# todo : norm spectra --> plot spectra --> print file


def spec_quick_view(ax, spec_chunks_list, wave_intervals=None,
                    wave_centers = None,
                    xtick_step = 50.,
                    num_spec_perpage=30,
                    offset_perspec=0.5,
                    wave_label='abs',
                    xtick_gap_fraction = .1,
                    cl='', *args, **kwargs):
    """
    Parameters
    ----------
    ax: axes
        the axes on which spectra will be plotted
    spec_chunks_list: list of list
        a list of spec chunks
    ranges: None or Nx2 array
        you have to specify ranges if you have more than one chunks for a spec

    """
    # for XTICK -----------------------------------------------------



    wave_intervals = np.array(wave_intervals)
    N_spec = len(spec_chunks_list)
    N_chunks = len(wave_intervals)

    wave_totalspan = np.sum(np.diff(wave_intervals))
    wave_gap = xtick_gap_fraction * wave_totalspan / (N_chunks - 1)
    print 'wave_gap', wave_gap

    # calculate how much wavelength should shift (to the left side)
    PLOT_START = 0.
    wave_offset_intercept = wave_intervals[0, 0] - PLOT_START
    wave_offset_spec = np.hstack((np.array(0), (wave_intervals[1:,0]-wave_intervals[:-1,1]).flatten()))
    wave_offset_gap = np.roll(np.ones(N_chunks) * wave_gap, 1)
    wave_offset_gap[0] = 0
    wave_offset = np.cumsum(wave_offset_spec - wave_offset_gap) + wave_offset_intercept
    wave_offset

    # for XTICK ------------------------------------------------------
    xtick_pos, xtick_lab = [], []
    for i_chunk in xrange(N_chunks):
        n_xtick_l = np.int(np.abs((wave_centers[i_chunk] - wave_intervals[i_chunk][0]) / xtick_step))
        n_xtick_r = np.int(np.abs((wave_centers[i_chunk] - wave_intervals[i_chunk][1]) / xtick_step))
        xtick_pos_, xtick_lab_ = generate_chunk_xtick_label(
            n_xtick_l, n_xtick_r, wave_centers[i_chunk],
            xtick_step=xtick_step, wave_offset=wave_offset[i_chunk])
        xtick_pos.extend(xtick_pos_)
        xtick_lab.extend(xtick_lab_)

    # PLOT -----------------------------------------------------------
    for i in xrange(N_spec):
        # for each list of spec chunks, do these:

        spec_chunks = spec_chunks_list[i]

        offset = -1 + (i+1) * offset_perspec

        for j in xrange(len(spec_chunks)):
            spec_chunk = spec_chunks[j]
            ax.plot(spec_chunk['wave']- wave_offset[j], spec_chunk['flux'] + offset, cl, *args, **kwargs)

    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(xtick_lab)
    return 0


def generate_chunk_xtick_label(n_xtick_l, n_xtick_r, wave_center,
                               xtick_step=50., wave_offset=0.):
    # relative xtick position
    xtick_pos_rel = np.arange(-n_xtick_l, n_xtick_r+1) * xtick_step

    # xtick labels
    xtick_lab = ['%.0f' % xtick_pos_rel_ for xtick_pos_rel_ in xtick_pos_rel]
    xtick_lab[n_xtick_l] = '%.1f' % wave_center

    # xtick positions
    xtick_pos = xtick_pos_rel + wave_center - wave_offset

    return xtick_pos, xtick_lab


def test_spec_quick_view():
    bc03dir = '/home/cham/PycharmProjects/bopy/bopy/data/model_bc03/'
    bc03cat = Table.read(
        '/home/cham/PycharmProjects/bopy/bopy/data/model_bc03/Base.BC03.N.csv')

    wave_intervals = [[5600, 5900], [6100., 6400]]

    specs = []
    for i in xrange(len(bc03cat)):
        fp = bc03dir + bc03cat['specfile'][i]
        spec = np.loadtxt(fp)
        spec = Spec(Table(spec, names=['wave', 'flux']))
        spec = norm_spec_pixel(spec, 5800.)
        specs.append(spec.extract_chunk_wave_interval(wave_intervals))

    fig = plt.figure('test spec_quick_view method')
    ax = fig.add_subplot(111)
    spec_quick_view(ax, specs, wave_intervals=wave_intervals,
                    wave_centers=[5780.1, 6284.1],
                    xtick_step = 50., xtick_gap_fraction = 0.02,
                    num_spec_perpage=30, offset_perspec=0.1, wave_label='abs',
                    cl='')
    # load a list of spectra

    # break into chunks

    # plot them

    print('@Cham: test spec_quick_view OK!')
    return 0


if __name__ == '__main__':
    test_spec_quick_view()

