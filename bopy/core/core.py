# -*- coding: utf-8 -*-
"""
@author: cham
Created on Mon Jan  4 20:08:23 2016
"""

import numpy as np



def logical_and_ext(ind_list):
    """ return the logical_and results for a list of index

    Parameters:
    -----------
    ind_list: list
        The list of ind arrays

    Returns:
    --------
    ind: array
        The logical_and result of the input ind_list
    """
    # assert that at least two ind are in the list
    assert len(ind_list) > 1
    # assert that all the ind have the same shape
    sp = ind_list[0].flatten().shape
    for ind in ind_list[1:]:
        assert sp == ind.flatten().shape
    # evaluate the logical_and result
    ind0 = ind_list[0].flatten()
    for ind in ind_list[1:]:
        ind0 = np.logical_and(ind0, ind.flatten())
    return ind0


def logical_or_ext(ind_list):
    """ return the logical_or results for a list of index

    Parameters:
    -----------
    ind_list: list
        The list of ind arrays
    
    Returns:
    --------
    ind: array
        The logical_or result of the input ind_list
    """
    # assert that at least two ind are in the list
    assert len(ind_list) > 1
    # assert that all the ind have the same shape
    sp = ind_list[0].flatten().shape
    for ind in ind_list[1:]:
        assert sp == ind.flatten().shape
    # evaluate the logical_and result
    ind0 = ind_list[0].flatten()
    for ind in ind_list[1:]:
        ind0 = np.logical_or(ind0, ind.flatten())
    return ind0

def core_print():
    print 'bopy is coming!'







if __name__ == '__main__':
    print logical_or_ext([np.ones(shape=(1,100)), np.zeros(shape=(2,50))])
    print logical_and_ext([np.ones(shape=(1,100)), np.zeros(shape=(2,50))])
    


