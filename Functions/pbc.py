import numpy as np


def pbc2d(r, box):
    '''
    :param r: ndarray of positions
    :param box: ndarray of box in cubic form: array([lx, ly, lz])
    :return: ndarray of positions
    '''
    return np.array([ x - box * np.round(x / box) for x in r ])
    #return r - box * np.round(r / box)

def pbc1d(r, box):
    return r - box*np.round(r/box)