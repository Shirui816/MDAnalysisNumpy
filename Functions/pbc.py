import numpy as np


def pbc(r, cbox):
    '''
    :param r: ndarray of positions
    :param box: ndarray of box in cubic form: array([lx, ly, lz])
    :return: ndarray of positions
    '''
    return np.array([ x - cbox * np.round(x / cbox) for x in r ])
    #return r - cbox * np.round(r / cbox)

def pbc1d(r, cbox):
    return r - cbox*np.round(r/cbox)