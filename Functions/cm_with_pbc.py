import numpy
from Functions.pbc import pbc2d, pbc1d
from cfunctions.cfunctions.functions import pbc2d as p2d
from cfunctions.cfunctions.functions import pbc as p1d


def cm(pos, box, mass=None):
    na = len(pos)
    if type(mass) == type(None):
        mass = numpy.ones(na)
        sum_mass = na
    else:
        sum_mass = sum(mass)
    pos_t = pbc2d(pos-pos[-1], box)
    pos_cm_t = sum((pos_t.T * mass).T)/sum_mass  #reduce to 1d
    pos_cm = pbc1d(pos_cm_t + pos[-1], box)
    return pos_cm

def cm_c(pos, box, mass=None):
    na = len(pos)
    if type(mass) == type(None):
        mass = numpy.ones(na)
        sum_mass = na
    else:
        sum_mass = sum(mass)
    pos_t = p2d(pos-pos[-1], box)
    pos_cm_t = sum((pos_t.T * mass).T)/sum_mass  #reduce to 1d
    pos_cm = pbc1d(pos_cm_t + pos[-1], box)
    return pos_cm
