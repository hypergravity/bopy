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
import astropy.constants as const
from astropy.table import Table
from bopy.spec.spec import spec_quick_init, Spec, norm_spec_pixel, wave2ranges

# input: a list of spectra
# todo : norm spectra --> plot spectra --> print file


def _spec_list_to_spec_chunk_list(n_spec, n_chunks, spec_list,
                                  norm_type=None, wave_intervals=None,
                                  q=0.90, delta_lambda=100.):

    # -----------------------------------------------------------------
    # norm_type:
    #   - None:             'flux' column directly used
    #   - 'continuum':      continuum normalization
    #   - 'median':         normalized to the median of whole spectrum
    #   - 'chunk_median:    normalized to the median of each chunk
    #   - 'pixel**':        normalized to a given pixel
    #                        useful for stellar population spectra
    # -----------------------------------------------------------------

    if norm_type is None or norm_type == 'chunk_median':
        # directly use the 'flux' column
        # or median normalization after breaking into chunks
        pass

    elif norm_type == 'continuum':  # PROBLEM!!!
        # 'continuum'
        for i in xrange(n_spec):
            spec_list[i].norm_spec_running_q(
                ranges=None,
                q=q, delta_lambda=delta_lambda, overwrite_flux=True)

    elif norm_type == 'median':
        # 'median'
        for i in xrange(n_spec):
            spec_list[i].norm_spec_median()

    elif norm_type[0:5] == 'pixel':
        # 'pixel'
        norm_wave = np.float64(norm_type[5:])
        for i in xrange(n_spec):
            spec_list[i].norm_spec_pixel(norm_wave)

    else:
        raise ValueError('norm_type wrong!')

    spec_chunk_list = [spec.extract_chunk_wave_interval(wave_intervals)
                       for spec in spec_list]

    if norm_type == 'chunk_median':
        # norm each spec chunk to its median
        for i in xrange(n_spec):
            for j in xrange(n_chunks):
                spec_chunk_list[i][j].norm_spec_median()

    return spec_chunk_list


def _calculate_wave_offset(n_spec, n_chunks, wave_intervals=None,
                           xtick_gap_fraction=0.05, plot_start=0.):
    assert wave_intervals is not None

    # total span of spec chunks
    wave_totalspan = np.sum(np.diff(wave_intervals))
    # gap is a fraction of total span
    wave_gap = xtick_gap_fraction * wave_totalspan / (n_chunks - 1)

    # calculate wave_offset, i.e.,
    # how much wavelength should shift (to the left side)
    wave_offset_intercept = wave_intervals[0, 0] - plot_start
    wave_offset_spec = np.hstack((np.array(0), (wave_intervals[1:,0]-wave_intervals[:-1,1]).flatten()))
    wave_offset_gap = np.roll(np.ones(n_chunks) * wave_gap, 1)
    wave_offset_gap[0] = 0
    wave_offset = np.cumsum(wave_offset_spec - wave_offset_gap) + wave_offset_intercept
    return wave_offset


def _xtick_pos_lab(n_chunks,
                   wave_centers, wave_intervals, wave_offset,
                   xtick_step,
                   xtick_format_str=('%.0f', '%.1f'),
                   xtick_label_type='wavelength'):
    # reshape wave_centers if necessary
    wave_centers = np.array(wave_centers)
    if wave_centers.ndim == 1:
        wave_centers = wave_centers.reshape((wave_centers.shape[0], 1))

    # initialize results
    xtick_pos, xtick_lab = [], []

    # xtick_label_type
    if xtick_label_type == 'velocity':
        # need to trasform 'wave_intervals' array to pseudo-volocity array
        wave_intervals = (wave_intervals - wave_centers) / wave_centers \
                         * const.c.value / 1000. + wave_centers

    for i_chunk in xrange(n_chunks):
        n_xtick_l = np.int(np.abs((wave_centers[i_chunk] - wave_intervals[i_chunk][0]) / xtick_step))
        n_xtick_r = np.int(np.abs((wave_centers[i_chunk] - wave_intervals[i_chunk][1]) / xtick_step))
        xtick_pos_, xtick_lab_ = _generate_chunk_xtick_label(
            n_xtick_l, n_xtick_r, wave_centers[i_chunk],
            xtick_step=xtick_step, wave_offset_chunk=wave_offset[i_chunk],
            xtick_format_str=xtick_format_str)
        xtick_pos.extend(xtick_pos_)
        xtick_lab.extend(xtick_lab_)
    return xtick_pos, xtick_lab


def _generate_chunk_xtick_label(n_xtick_l, n_xtick_r,
                                wave_center,
                                xtick_step=50.,
                                wave_offset_chunk=0.,
                                xtick_format_str=('%.0f', '%.1f')):
    # relative xtick position
    xtick_pos_rel = np.arange(-n_xtick_l, n_xtick_r+1) * xtick_step

    # xtick labels
    xtick_lab = [xtick_format_str[0] % xtick_pos_rel_ for xtick_pos_rel_ in xtick_pos_rel]
    xtick_lab[n_xtick_l] = xtick_format_str[1] % wave_center

    # xtick positions
    xtick_pos = xtick_pos_rel + wave_center - wave_offset_chunk

    return xtick_pos, xtick_lab


def spec_quick_view(ax,
                    spec_list,
                    norm_type='median',
                    q=0.90,
                    delta_lambda=100.,
                    wave_intervals=None,
                    wave_centers = None,
                    xtick_modify=True,
                    xtick_step = 50.,
                    xtick_gap_fraction = .1,
                    xtick_start=0.,
                    xtick_format_str=('%.0f', '%.1f'),
                    xtick_label_type='wavelength',
                    offset_perspec=0.5,
                    wave_label='abs',
                    cl='',
                    *args, **kwargs):
    """
    Parameters
    ----------
    ax: axes
        the axes on which spectra will be plotted
    spec_list: list of list
        a list of spec
    ranges: None or Nx2 array
        you have to specify ranges if you have more than one chunks for a spec

    """

    wave_intervals = np.array(wave_intervals)
    n_chunks = len(wave_intervals)
    n_spec = len(spec_list)

    # norm spec & break spec into chunks -------------------------------------

    # norm_type:
    #   - None:             'flux' column directly used
    #   - 'continuum':      continuum normalization
    #   - 'median':         normalized to the median of whole spectrum
    #   - 'chunk_median:    normalized to the median of each chunk
    #   - 'pixel**':        normalized to a given pixel
    #                        useful for stellar population spectra

    spec_chunk_list = _spec_list_to_spec_chunk_list(
        n_spec, n_chunks, spec_list, norm_type, wave_intervals,
        q=q, delta_lambda=delta_lambda)

    # for XTICK --------------------------------------------------------------

    offset_specs = np.arange(n_spec) * offset_perspec

    if not xtick_modify:
        # just plot spec chunks, don't modify XTICK & XTICKLABELS
        for i in xrange(n_spec):
            # for each list of spec chunks, do these:
            spec_chunks = spec_chunk_list[i]
            for j in xrange(n_chunks):
                spec_chunk = spec_chunks[j]
                ax.plot(spec_chunk['wave'],
                        spec_chunk['flux'] + offset_specs[i],
                        cl, *args, **kwargs)

        # set xlim
        ax.set_xlim(np.min(wave_intervals[:, 0]), np.max(wave_intervals[:, 1]))
    else:
        assert wave_intervals is not None
        # need to modify XTICK & XTICKLABELS

        # 1> calculate wave_offset
        wave_offset = _calculate_wave_offset(n_spec, n_chunks, wave_intervals,
                                             xtick_gap_fraction, xtick_start)
        print 'wave_offset', wave_offset
        for i in xrange(n_spec):
            # for each list of spec chunks, do these:
            spec_chunks = spec_chunk_list[i]
            for j in xrange(n_chunks):
                spec_chunk = spec_chunks[j]
                ax.plot(spec_chunk['wave'] - wave_offset[j],  # wave_offset
                        spec_chunk['flux'] + offset_specs[i],
                        cl, *args, **kwargs)

        # 2> calculate xtick position & labels
        xtick_pos, xtick_lab = _xtick_pos_lab(
            n_chunks, wave_centers, wave_intervals, wave_offset, xtick_step,
            xtick_format_str=xtick_format_str,
            xtick_label_type=xtick_label_type)

        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(xtick_lab)

        # set xlim
        ax.set_xlim(
            np.min(wave_intervals[:, 0].flatten()-wave_offset.flatten()),
            np.max(wave_intervals[:, 1].flatten()-wave_offset.flatten()))
    return ax, spec_chunk_list, wave_offset, offset_specs


def test_spec_quick_view():
    # read BC03 catalog
    print('-------------------------------------------------------')
    print('@Cham: reading BC03 catalog ...')
    bc03cat = Table.read(
        '/home/cham/PycharmProjects/bopy/bopy/data/model_bc03/Base.BC03.N.csv')
    bc03dir = '/home/cham/PycharmProjects/bopy/bopy/data/model_bc03/'

    # prepare wave_intervals
    wave_intervals = [[5550, 5800], [6250., 6300]]

    # break spectra into chunks
    print('@Cham: loading BC03 spectra and break into chunks ...')
    spec_list = []
    for i in xrange(10):
        fp = bc03dir + bc03cat['specfile'][i]
        spec = np.loadtxt(fp)
        spec = Spec(Table(spec, names=['wave', 'flux']))
        # spec = norm_spec_pixel(spec, 6000.)
        # spec_chunks_list.append(
        #     spec.extract_chunk_wave_interval(wave_intervals))
        spec_list.append(spec)
    # plot spec chunks
    print('@Cham: ploting spec chunks ...')
    fig = plt.figure('test spec_quick_view method', figsize=(10, 20))
    ax = fig.add_subplot(111)
    ax, spec_chunk_list, wave_offset, offset_specs = \
        spec_quick_view(ax, spec_list,
                        norm_type='chunk_median',
                        wave_intervals=wave_intervals,
                        wave_centers=[5780.1, 6284.1],
                        xtick_step=50.,
                        xtick_gap_fraction=0.02,
                        xtick_label_type='velocity',
                        # num_spec_perpage=30,
                        offset_perspec=0.5,
                        wave_label='abs',
                        cl='k-')
    print spec_chunk_list[0][0]
    print('@Cham: saving figure ...')
    fig.savefig('/home/cham/PycharmProjects/bopy/bopy/data/test_spec_quick_view/test.pdf')
    print('@Cham: test spec_quick_view OK!')
    return 0


if __name__ == '__main__':
    test_spec_quick_view()

