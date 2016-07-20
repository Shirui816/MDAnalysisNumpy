class position(object):
    def __cinit__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


from libc.math cimport round, sqrt, floor

def pbc(double r, double d):
    return r - d * round(r/d)


def distance(double x0, double y0, double z0, double x1, double y1, double z1, double Lx, double Ly, double Lz):
    return sqrt(pbc(x1 - x0, Lx) ** 2 + pbc(y1 - y0, Ly) ** 2 + pbc(z1 - z0, Lz) ** 2)


def icell(long ix, long iy, long iz, long ibx, long iby, long ibz):
    return (ix + ibx) % ibx + ((iy + iby) % iby) * ibx + ((iz + ibz) % ibz) * iby * ibx


def wcell(double x, double y, double z, ib, box):
    return int((x/box.x + 0.5) * ib.x) + int((y/box.y + 0.5) * ib.y) * ib.x + int((z/box.z + 0.5) * ib.z) * ib.y * ib.x

import numpy as np
cimport cython
cimport numpy as np

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def pbc2d(np.ndarray[double, ndim=2] pos, np.ndarray[double, ndim=1] box):
    cdef long na, d, i
    na = pos.shape[0]
    d = pos.shape[1]
    cdef np.ndarray res = np.zeros([na, d], dtype=np.float)
    for i in range(na):
        res[i, 0] = pbc(pos[i, 0], box[0])
        res[i, 1] = pbc(pos[i, 1], box[1])
        res[i, 2] = pbc(pos[i, 2], box[2])
    return res

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def cm_cc(np.ndarray[double, ndim=2] pos, np.ndarray[double, ndim=1] box):
    cdef long na, d, i
    cdef double xcm=0
    cdef double ycm =0
    cdef double zcm=0
    cdef np.ndarray res = np.zeros(3, dtype=np.float)
    na = pos.shape[0]
    d = pos.shape[1]
    cdef np.ndarray dp = np.zeros([na, d], dtype=np.float)
    for i in range(na-1):
        dp[i+1, 0] = pbc(pos[i+1,0] - pos[i, 0], box[0]) + dp[i,0]
        dp[i+1, 1] = pbc(pos[i+1,1] - pos[i, 1], box[1]) + dp[i,1]
        dp[i+1, 2] = pbc(pos[i+1,2] - pos[i, 2], box[2]) + dp[i,2]
    for i in range(na):
        xcm += dp[i,0]
        ycm += dp[i,1]
        zcm += dp[i,2]
    res[0] = pbc(xcm/na + pos[0,0], box[0])
    res[1] = pbc(ycm/na + pos[0,1], box[1])
    res[2] = pbc(zcm/na + pos[0,2], box[2])
    return res




#@cython.boundscheck(False) # turn off bounds-checking for entire function
#@cython.wraparound(False)  # turn off negative index wrapping for entire function
def RgRadial(np.ndarray r_pos, np.ndarray box):
        cdef long m = r_pos.shape[0]
        cdef double rg2 = 0
        cdef double rgn2 = 0
        cdef long i, j
        for i in range(m-1):
                for j in range(i+1, m):
                        #d = pbc1d((r_pos[i]+r_pos[j])/2, box)
                        dx = pbc((r_pos[i,0] + r_pos[j,0])/2, box[0])
                        dy = pbc((r_pos[i,1] + r_pos[j,1])/2, box[1])
                        dz = pbc((r_pos[i,2] + r_pos[j,2])/2, box[2])
                        sx = pbc(r_pos[i][0]-r_pos[j][0], box[0])
                        sy = pbc(r_pos[i][1]-r_pos[j][1], box[1])
                        sz = pbc(r_pos[i][2]-r_pos[j][2], box[2])
                        #s = pbc1d(r_pos[i]-r_pos[j], box)
                        dr = dx*dx+dy*dy+dz*dz
                        rg2 += sx*sx+sy*sy+sz*sz
                        rgn2 += (sx*dx+sy*dy+sz*dz)**2 / dr
        return(rg2/m**2, rgn2/m**2)



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def RgRadial_eff_idx(np.ndarray[double, ndim=2] r_pos, np.ndarray[double, ndim=1] box):
        cdef long m = r_pos.shape[0]
        cdef double rg2 = 0
        cdef double rgn2 = 0
        cdef long i
        cdef long j
        cdef double dx, dy, dz, sx, sy, sz, dr
        for i in range(m-1):
                for j in range(i+1, m):
                        #d = pbc1d((r_pos[i]+r_pos[j])/2, box)
                        dx = pbc((r_pos[i,0] + r_pos[j,0])/2, box[0])
                        dy = pbc((r_pos[i,1] + r_pos[j,1])/2, box[1])
                        dz = pbc((r_pos[i,2] + r_pos[j,2])/2, box[2])
                        sx = pbc(r_pos[i,0]-r_pos[j,0], box[0])
                        sy = pbc(r_pos[i,1]-r_pos[j,1], box[1])
                        sz = pbc(r_pos[i,2]-r_pos[j,2], box[2])
                        #s = pbc1d(r_pos[i]-r_pos[j], box)
                        dr = dx*dx+dy*dy+dz*dz
                        rg2 += sx*sx+sy*sy+sz*sz
                        rgn2 += (sx*dx+sy*dy+sz*dz)**2 / dr
        return(rg2/m**2, rgn2/m**2)
