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
- Wed Mar  30 12:17:00 2016     output random dots on an image

Modifications
-------------
- Wed Mar  30 12:17:00 2016     output random dots on an image

Aims
----
- Wed Mar  30 12:17:00 2016     output random dots on an image

"""

import numpy as np


def random_dots(immask, n_dots=10, mindist=None, n_loop_max=1E9):
    """
    Parameters
    ----------
    #imdata: numpy.ndarray (float)
    #    image data
    imdata: numpy.ndarray (bool)
        image mask data.
    n_dots: int
        number of requested random dots
    mindist: None|tuple of float
        if None, no criteria on distance between dots
        if specified, the min distance between dots is set
    n_loop_max: int
        max loop can be used to generate random dots

    Returns
    -------
    x, y: X and Y index

    """
    assert np.ndim(immask) == 2
    assert n_dots > 0

    if mindist is None:
        # no criteria on min distance
        sz = immask.shape
        n_ele = sz[0] * sz[1]
        n_val = np.float(np.sum(immask))
        n_est = np.int64(10 * n_dots / (n_val / n_ele))

        x = np.random.randint(0, sz[1], n_est)
        y = np.random.randint(0, sz[0], n_est)
        ind_val = immask[y, x]
        sum_val = np.sum(ind_val)
        if sum_val > n_dots:
            # successful
            return x[ind_val][:n_dots], y[ind_val][:n_dots], True
        else:
            # unsucessful
            return random_dots(immask, n_dots=n_dots, mindist=mindist)
    else:
        # set min distance
        assert len(mindist) == 2
        mindist_y, mindist_x = mindist
        sz = immask.shape

        x, y = np.zeros((n_dots,), np.int64), np.zeros((n_dots,), np.int64)
        x[0], y[0] = _random_dot_valid(immask, sz)
        c = 1
        n_loop = 0
        print('@Cham: got random dots [%s/%s] (y, x) = (%d,%d),  ...' % (c, n_dots, y[0], x[0]))
        while c < n_dots and n_loop < n_loop_max:
            n_loop += 1
            x_, y_ = _random_dot_valid(immask, sz)
            dist = ((x_-x[:c])/mindist_x)**2. + ((y_-y[:c])/mindist_y)**2.
            if np.all(dist.flatten() > 1.):
                x[c], y[c] = x_, y_
                c += 1
                print('@Cham: got random dots [%s/%s] (y, x) = (%d,%d),  ...' % (c, n_dots, y_, x_))
        if c < n_dots:
            print('@Cham: failed to get all requested random dots!')
            return x, y, False
        else:
            print('@Cham: got all requested random dots successfully!')
            return x, y, True


def _random_dot_valid(immask, sz=None):
    assert np.sum(immask) > 0
    if sz is None:
        sz = immask.shape
    x_, y_ = np.random.randint(0, sz[1]), np.random.randint(0, sz[0])
    while not immask[y_, x_]:
        x_, y_ = np.random.randint(0, sz[1]), np.random.randint(0, sz[0])
    return x_, y_


def _test():
    imdata = np.array([[1,2,3],[4,5,6]])
    immask = imdata>2
    print imdata
    print immask
    print imdata[immask]
    print np.ndim(immask)


def _test2():
    imdata = np.ones((700,700))
    for i in xrange(700):
        for j in xrange(700):
            if i < 200:
                imdata[i, j] = np.nan
    immask = np.logical_not(np.isnan(imdata))
    print np.sum(immask)
    x,y,flag = random_dots(immask, 30, (5., 5.),1E9)
    print x,y,flag


if __name__ == '__main__':
    _test2()