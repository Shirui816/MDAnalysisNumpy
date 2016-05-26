import numpy
from Functions.cm_with_pbc import cm
from Functions.pbc import pbc2d, pbc1d


def shells(center, pos, box, s = 300):
        '''
        Centeroid position is needed; 300 shells for default
        '''
        n, m = pos.shape # Ensures that pos is a 2D array
        result = numpy.zeros((n, s), dtype=numpy.int)
        L = box.min()/2
        delta = L/s
        pbc_pos = pbc2d(pos - center, box)
        for i in range(n):
                d = sum(pbc_pos[i]**2)**0.5
                if not d > L:
                        result[i][int(d/delta)] += 1
        return(result)

def shells_c0(pos_c0, box, s = 300):
        '''
        Positions already centered at 0; 300 shells for default
        '''
        n, m = pos_c0.shape # Ensures that pos is a 2D array
        result = numpy.zeros((n, s), dtype=numpy.int) # N paritcles with s shells
        L = box.min()/2
        delta = L/s
        for i in range(n):
                d = sum(pos_c0[i]**2)**0.5
                if not d > L:
                        result[i][int(d/delta)] += 1
        return(result)