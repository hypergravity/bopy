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
- Fri Mar  25 17:57:00 2016     get a set of PADOVA isochrones grid

Modifications
-------------
- Fri Mar  25 17:57:00 2016     get a set of PADOVA isochrones grid

Aims
----
- get a set of PADOVA isochrones grid
- output the combined isochrone table

"""

import numpy as np
from bopy.helpers.ezpadova import cmd
from astropy.table import Table, vstack, Column
from astropy.io import fits
from scipy.interpolate import PchipInterpolator
from joblib import Parallel, delayed


Zsun = 0.0152  # this value from CMD website
Zmin = 0.0001
Zmax = 0.07
logtmax = 10.13
logtmin = 1.


def _find_valid_grid(grid_feh, grid_logt, Zsun=0.0152):
    """ return valid grid of [feh, logt] """
    grid_feh = np.array(grid_feh).flatten()
    grid_logt = np.array(grid_logt).flatten()
    grid_Z = 10. ** grid_feh * Zsun
    ind_valid_Z = np.logical_and(grid_Z >= Zmin, grid_Z <= Zmax)
    ind_valid_logt = np.logical_and(grid_logt >= logtmin, grid_logt <= logtmax)
    vgrid_feh = grid_feh[ind_valid_Z]
    vgrid_logt = grid_logt[ind_valid_logt]
    print('@Cham: -----------------------------------------------------------')
    print('@Cham: the valid range for Z & logt are (%s, %s) & (%s, %s).'
          % (Zmin, Zmax, logtmin, logtmax))
    print('@Cham: -----------------------------------------------------------')
    print('@Cham: valid input feh are:')
    print(grid_feh[ind_valid_Z])
    print('@Cham: valid input logt are:')
    print(grid_logt[ind_valid_logt])
    print('@Cham: INvalid input feh are:')
    print(grid_feh[~ind_valid_Z])
    print('@Cham: INvalid input logt are:')
    print(grid_logt[~ind_valid_logt])
    print('@Cham: -----------------------------------------------------------')
    return vgrid_feh, vgrid_logt


def get_isochrone_grid(grid_feh,
                       grid_logt,
                       model='parsec12s',
                       phot='sloan',
                       Zsun=0.0152,
                       parflag=True,
                       n_jobs=8,
                       **kwargs):
    """ get a list of isochrones using EZPADOVA

    Parameters
    ----------
    grid_feh: array
        [Fe/H] grid
    grid_logt: array
        logt grid
    model: string
        default is 'parsec12s'
    phot: string
        default is 'sloan'
    Zsun: float
        default is 0.0152
    parflag: bool
        default is True
        if True, use JOBLIB to get isochrones in parallel
    n_jobs: int
        if parflat is True, specify number of jobs in JOBLIB

    Returns
    -------
    vgrid_feh, vgrid_logt, isoc_list, grid_list

    """
    # validate grid
    vgrid_feh, vgrid_logt = _find_valid_grid(grid_feh, grid_logt, Zsun=Zsun)

    # construct list
    grid_list = []
    for grid_feh_ in vgrid_feh:
        for grid_logt_ in vgrid_logt:
            grid_list.append((10.**grid_logt_, 10.**grid_feh_*Zsun))

    print('@Cham: you have requested for %s isochrones!')
    print('@Cham: -----------------------------------------------------------')

    # get isochrones
    if parflag:
        # get isochrones in parallel
        isoc_list = Parallel(n_jobs=n_jobs, verbose=True)(delayed(cmd.get_one_isochrone)(
            grid_list_[0], grid_list_[1], model=model, phot=phot, **kwargs) for grid_list_ in grid_list)
    else:
        # get isochrones sequentially
        isoc_list = []
        for i in xrange(len(grid_list)):
            grid_list_ = grid_list[i]
            print('@Cham: sending request for isochrone (logt=%s, [Fe/H]=%s) (t=%s, Z=%s) [%s/%s]...'
                  % (np.log10(grid_list_[0]), np.log10(grid_list_[1]/Zsun), grid_list_[0], grid_list_[1], i+1, len(grid_list)))
            isoc_list.append(
                Table(cmd.get_one_isochrone(grid_list_[0], grid_list_[1], model=model, phot=phot, **kwargs).data))
    print('@Cham: got all requested isochrones!')
    print('@Cham: -----------------------------------------------------------')
    print('@Cham: colnames are:')
    print(isoc_list[0].colnames)
    print('@Cham: -----------------------------------------------------------')
    return vgrid_feh, vgrid_logt, isoc_list, grid_list


def interpolate_to_cube(grid_feh, grid_logt, isoc_list, grid_list, grid_mini, cube_quantities=[]):
    """ interpolate a slit of isochrones into data cubes
    grid_feh: array
        [Fe/H] grid
    grid_logt: array
        logt grid
    isoc_list: list
        list of isochrones (in astropy.table.Table form)
    grid_list: list
        a list of (logt, Z) tuples corresponding to isoc_list
    grid_mini: array
        the M_ini array to which interpolate into
    cube_quantities: list
        a list of names of the quantities to be interpolated

    Returns
    -------
    cube_data_list, cube_name_list

    """
    # flatten grid
    grid_feh = grid_feh.flatten()
    grid_logt = grid_logt.flatten()
    grid_mini = grid_mini.flatten()

    # mesh cube
    cube_logt, cube_feh, cube_mini, = np.meshgrid(grid_logt, grid_feh, grid_mini)
    cube_size = cube_feh.shape
    print('@Cham: cube shape: ', cube_size)
    print('@Cham: -----------------------------------------------------------')

    # determine cube-quantities
    if len(cube_quantities) == 0:
        # all the quantities besides [feh, logt, Mini]
        colnames = list(isoc_list[0].colnames)
        assert colnames[0] == 'Z'
        assert colnames[1] == 'logageyr'
        assert colnames[2] == 'M_ini'
        cube_quantities = colnames[3:]
        print('@Cham: I will now begin interpolating these quantities into cubes ...')
        print(cube_quantities)
        print('@Cham: -----------------------------------------------------------')
    else:
        print('@Cham: I will now begin interpolating these quantities into cubes ...')
        print(cube_quantities)
        print('@Cham: -----------------------------------------------------------')

    # smoothing along M_ini
    for i in xrange(len(isoc_list)):
        # Tablize
        if not isinstance(isoc_list[i], Table):
            isoc_list[i] = Table(isoc_list[i].data)
        print '@Cham: smoothing isochrones [%s/%s] ...' % (i+1, len(isoc_list))
        # smoothing M_ini
        ind_same_mini = np.hstack((False, np.diff(isoc_list[i]['M_ini'].data)==0))
        sub_same_mini = np.arange(len(isoc_list[i]))[ind_same_mini]
        isoc_list[i].remove_rows(sub_same_mini)
        print '@Cham: removing %s rows for this table ...' % len(sub_same_mini)
    print('@Cham: -----------------------------------------------------------')

    # interpolation
    cube_data_list = [cube_feh, cube_logt, cube_mini]
    cube_name_list = ['feh', 'logt', 'M_ini']
    for k in xrange(len(cube_quantities)):
        cube_name = cube_quantities[k]
        c = 0
        cube_data = np.ones(cube_size) * np.nan
        for i in xrange(len(grid_feh)):
            for j in xrange(len(grid_logt)):
                this_isoc = isoc_list[c]
                P = PchipInterpolator(this_isoc['M_ini'].data, this_isoc[cube_name].data, extrapolate=False)
                cube_data[i, j, :] = P(grid_mini)
                print('@Cham: interpolating cube quantity [%s] {quantity: %s/%s} (%s/%s) ...'
                      % (cube_name, k+1, len(cube_quantities), c+1, len(grid_feh)*len(grid_logt)))
                c += 1
        cube_data_list.append(cube_data)
        cube_name_list.append(cube_name)
    print('@Cham: -----------------------------------------------------------')

    return cube_data_list, cube_name_list


def cubelist_to_hdulist(cube_data_list, cube_name_list):
    """ transform data cubes into fits HDU list """
    print('@Cham: transforming data cubes into HDU list ...')
    # construct Primary header
    header = fits.Header()
    header['author'] = 'Bo Zhang (@NAOC)'
    header['data'] = 'isochrone cube'
    header['software'] = 'cube constructed using BOPY'

    # initialize HDU list
    hl = [fits.hdu.PrimaryHDU(header=header)]

    # construct HDU list
    for i in xrange(len(cube_data_list)):
        hl.append(fits.hdu.ImageHDU(data=cube_data_list[i], name=cube_name_list[i]))

    print('@Cham: -----------------------------------------------------------')
    return fits.HDUList(hl)


def combine_isochrones(isoc_list):
    """ combine isochrone Tables into 1 Table"""
    if isinstance(isoc_list[0], Table):
        # assume that these data are all Table
        comb_isoc = vstack(isoc_list)
    else:
        for i in xrange(isoc_list):
            isoc_list[i] = Table(isoc_list[i])
        comb_isoc = vstack(isoc_list)
    return comb_isoc


def write_isoc_list(isoc_list,
                    grid_list,
                    dirpath='comb_isoc_parsec12s_sloan',
                    extname='.fits',
                    Zsun=0.0152):
    """ write isochrone list into separate tables """
    assert len(isoc_list) == len(grid_list)
    for i in xrange(len(isoc_list)):
        fp = dirpath + \
            '_ZSUN' + ('%.5f' % Zsun).zfill(7) + \
            '_LOGT' + ('%.3f' % np.log10(grid_list[i][0])).zfill(6) + \
            '_FEH' + ('%.3f' % np.log10(grid_list[i][1]/Zsun)).zfill(6) + \
            extname
        print('@Cham: writing table [%s] [%s/%s]...' % (fp, i+1, len(isoc_list)))
        isoc_list[i].write(fp, overwrite=True)
    return


def _test():
    # set grid
    grid_logt = [6, 7., 9]
    grid_feh  = [-2.2, -1., 0, 1., 10]
    grid_mini = np.arange(0.01, 12, 0.01)

    # get isochrones
    vgrid_feh, vgrid_logt, isoc_list, grid_list = get_isochrone_grid(
        grid_feh, grid_logt, model='parsec12s', phot='sloan', parflag=True)

    # transform into cube data
    cube_data_list, cube_name_list = interpolate_to_cube(
        vgrid_feh, vgrid_logt, isoc_list, grid_list, grid_mini,
        cube_quantities=['M_act', 'g', 'r'])

    # cube HDUs
    hl = cubelist_to_hdulist(cube_data_list, cube_name_list)
    hl.info()
    # hl.writeto()

    # combine isochrone tables
    comb_isoc = combine_isochrones(isoc_list)
    # comb_isoc.write()

    # write isochrone list into separate files
    # write_isoc_list(isoc_list, grid_list, '/pool/comb_isoc')
    return hl


if __name__ == '__main__':
    _test()






















