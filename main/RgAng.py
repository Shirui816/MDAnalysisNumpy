#!/bin/python
from Functions.cm_with_pbc import cm
from Functions.pbc import pbc2d, pbc1d
import numpy
from MoleculeClassify.hoomd_mols import hoomd_mols
from math import acos, pi, sqrt


from cfunctions.cfunctions.functions import pbc, RgRadial

def RgTensor(pos, box):
        n, m = pos.shape
        assert(m==3)
        center = cm(pos, box)
        pos = pbc2d(pos - center, box)
        x = pos.T[0]
        y = pos.T[1]
        z = pos.T[2]
        Sxx = sum(x*x)/n
        Syy = sum(y*y)/n
        Szz = sum(z*z)/n
        Sxy = Syx = sum(x*y)/n
        Sxz = Szx = sum(x*z)/n
        Syz = Szy = sum(y*z)/n
        return numpy.array([[Sxx, Sxy, Sxz],[Syx, Syy, Syz],[Szx, Szy, Szz]])

def MajorAxis(rg_tensor):
        w, v = numpy.linalg.eig(rg_tensor)
        vl = [ x for y, x in sorted(zip(w, v.T), reverse=True) ]
        return vl[0]


def MajorAxis_np(rg_tensor):
        w, v = numpy.linalg.eig(rg_tensor)
        inds = w.argsort()
        vl = v.T[inds]
        return vl[-1]
def RgRadial_py(r_pos, box):
        m, n = r_pos.shape
        rg2 = rgn2 = 0
        for i in range(m-1):
                for j in range(i, m):
                        #d = pbc1d((r_pos[i]+r_pos[j])/2, box)
                        dx = pbc(r_pos[i][0] + r_pos[j][0]/2, box[0])
                        dy = pbc(r_pos[i][1] + r_pos[j][1]/2, box[1])
                        dz = pbc(r_pos[i][2] + r_pos[j][2]/2, box[2])
                        d = np.array([dx, dy,dz])
                        sx = pbc(r_pos[i][0]-r_pos[j][0], box[0])
                        sy = pbc(r_pos[i][1]-r_pos[j][1], box[1])
                        sz = pbc(r_pos[i][2]-r_pos[j][2], box[2])
                        #s = pbc1d(r_pos[i]-r_pos[j], box)
                        s = np.array([sx, sy, sz])
                        dr = d.dot(d)
                        rg2 += s.dot(s)
                        rgn2 += s.dot(d) ** 2 / dr
        return(rg2, rgn2)

def shell_id(r_pos, binsize):
        sid = 0
        m, n = r_pos.shape
        for p in r_pos:
                r = sqrt(p.dot(p))
                sid += int(r/binsize)
        sid /= m
        return(sid)

def main(seg, xml, binsize = 0.1):
        pos = xml.nodes['position']
        box = xml.box
        bins = int(0.5 * sqrt(box.dot(box))/binsize) + 1
        rg2s = numpy.zeros((bins,), dtype=numpy.float)
        rgn2s = numpy.zeros((bins,),dtype=numpy.float)
        angs = numpy.zeros((bins,),dtype=numpy.float)
        cid = numpy.zeros((bins,),dtype=numpy.float)
        types = xml.nodes['type']
        posA = pos[types==b'A']
        cmA = cm(posA, xml.box)
        r_pos = pbc2d(pos-cmA, box)
        mols = hoomd_mols(xml)
        countSeg = 0
        for mol in mols.mol_idxes[1:]:
                length = len(mol)
                for s in range(0, length, seg):
                        if s + seg > length:
                                continue
                        countSeg += 1
                        seg_idx = mol[s:s+seg]
                        seg_pos = r_pos[seg_idx]
                        cm_seg = cm(seg_pos, box)
                        r_cm = sqrt(cm_seg.dot(cm_seg))
                        sid = shell_id(seg_pos, binsize)
                        cid[sid] += 1
                        rg2, rgn2 = RgRadial(seg_pos, box)
                        rgtensor = RgTensor(seg_pos, box)
                        maxis = MajorAxis_np(rgtensor)
                        r_maixs = sqrt(maxis.dot(maxis))
                        angle = maxis.dot(cm_seg)/r_cm/r_maixs
                        angs[sid] += acos(angle)
                        rg2s[sid] += rg2
                        rgn2s[sid] += rgn2
        for i in range(bins):
                if cid[i] ==0:
                        cid[i] = 1
        rg2s/=cid
        rgn2s/=cid
        angs/=cid
        angs = angs/pi * 180
        rgfile = open('rg.txt', 'w')
        rgfile.write('#r, rg2, rgn2\n')
        angfile = open('ang.txt', 'w')
        angfile.write('#r, ang\n')
        i = 0
        for rg2, rgn2 in zip(rg2s, rgn2s):
                r = binsize * (0.5+ i)
                rgfile.write("%.4f %.4f %.4f\n" % (r, rg2, rgn2))
                i += 1
        for i, ang in enumerate(angs):
                angfile.write("%.4f %.4f\n" % (binsize * (i+0.5), ang))
        rgfile.close()
        angfile.close()