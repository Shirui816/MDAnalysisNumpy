import numpy as np

def cm(pos, mass=None):
    na = len(pos)
    if type(mass)==type(None):
        mass = np.ones(na)
        sum_mass = na
    else:
        sum_mass = sum(mass)
    return sum((pos.T*mass).T)/sum_mass