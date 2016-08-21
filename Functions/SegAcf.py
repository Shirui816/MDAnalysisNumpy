#!/bin/python
from Functions.cm_with_pbc import cm, cm_c
from Functions.pbc import pbc1d
import numpy
from MoleculeClassify.hoomd_mols import hoomd_mols
from math import acos, pi, sqrt
from cfunctions.cfunctions.functions import cm_cc, pbc2d

from cfunctions.cfunctions.functions import pbc, RgRadial, RgRadial_eff_idx


def shell_id(r_pos, binsize):
	sid = 0
	m, n = r_pos.shape
	for p in r_pos:
		r = sqrt(p.dot(p))
		sid += int(r/binsize)
	sid /= m
	return(int(sid))

def shell_dist(r_pos, binsize, bins):
	sdist = numpy.zeros((bins,), dtype=numpy.float)
	m, n = r_pos.shape
	for p in r_pos:
		r = sqrt(p.dot(p))
		sdist[int(r/binsize)] += 1
	return(sdist, sdist/sum(sdist))

def shell_cm(r_pos, binsize, box):
	cmv = cm_cc(r_pos, box)
	cmr = sqrt(cmv.dot(cmv))
	return(int(cmr/binsize))

def main(segs, xml, mols, binsize = 0.1, fn='', COMP=1):
	pos = xml.nodes['position']
	box = xml.box
	bins = int(0.5 * sqrt(box.dot(box))/binsize) + 1
	start = len(mols.sbody)
	if not COMP:
		bc = mols.bcs()[0]
	#mols = hoomd_mols(xml)
	r_pos=pos
	if not COMP:
        	r_pos = pbc2d(pos - bc, box)
	countSeg = {}
	fs = {}
	for seg in segs:
		countSeg[seg] = 0
	shells = numpy.zeros((bins,))
	for mol in mols.mol_idxes[start:]:
		length = len(mol)
		for seg in segs:
			fs[seg] = open('%s_ree.txt' % (seg), 'a')
			for s in range(0, length, seg):
				if s + seg > length:
					continue
				countSeg[seg] += 1
				seg_idx = mol[s:s+seg]
				seg_pos = r_pos[seg_idx]
				sid = 0
				if not COMP:
					sid = shell_cm(seg_pos, binsize, box)
				reevx = pbc(seg_pos[-1][0] - seg_pos[0][0], box[0])
				reevy = pbc(seg_pos[-1][1] - seg_pos[0][1], box[1])
				reevz = pbc(seg_pos[-1][2] - seg_pos[0][2], box[2])
				#ss, sdist = shell_dist(seg_pos, binsize, bins)
				#shells += s
				#sid = shell_id(seg_pos, binsize)
				fs[seg].write("%.4f %.4f %.4f %.4f\n" % (reevx, reevy, reevz, sid))
			fs[seg].close()
	#for i in range(bins):
		#if cid[i] ==0:
			#cid[i] = 1
	return(countSeg)
	
	#rgfile = open(fn+'rg.txt', 'w')
	#rgfile.write('#r, rg2, rgn2\n')
	#angfile = open(fn+'ang.txt', 'w')
	#angfile.write('#r, ang\n')
	#i = 0
	#for rg2, rgn2 in zip(rg2s, rgn2s):
	#	r = binsize * (0.5+ i)
	#	rgfile.write("%.4f %.4f %.4f\n" % (r, rg2, rgn2))
	#	i += 1
	#for i, ang in enumerate(angs):
	#	angfile.write("%.4f %.4f\n" % (binsize * (i+0.5), ang))
	#rgfile.close()
	#angfile.close()





def main4pp(seg, xml, mols, binsize = 0.1, fn=''):
	pos = xml.nodes['position']
	box = xml.box
	bins = int(0.5 * sqrt(box.dot(box))/binsize) + 1
	rg2s = numpy.zeros((bins,), dtype=numpy.float)
	rgn2s = numpy.zeros((bins,),dtype=numpy.float)
	angs = numpy.zeros((bins,),dtype=numpy.float)
	cid = numpy.zeros((bins,),dtype=numpy.float)
	types = xml.nodes['type']
	posA = pos[types=='A']
	cmA = cm_cc(posA, xml.box)
	r_pos = pbc2d(pos-cmA, box)
	#mols = hoomd_mols(xml)
	countSeg = 0
	res = []
	shells = numpy.zeros((bins,))
	for mol in mols.mol_idxes[1:]:
	        length = len(mol)
	        for s in range(0, length, seg):
	                if s + seg > length:
	                        continue
	                countSeg += 1
	                seg_idx = mol[s:s+seg]
	                seg_pos = r_pos[seg_idx]
	                sid = shell_cm(seg_pos, binsize, box)
	                reevx = pbc(seg_pos[-1][0] - seg_pos[0][0], box[0])
	                reevy = pbc(seg_pos[-1][1] - seg_pos[0][1], box[1])
	                reevz = pbc(seg_pos[-1][2] - seg_pos[0][2], box[2])
	                #ss, sdist = shell_dist(seg_pos, binsize, bins)
	                #shells += s
	                #sid = shell_id(seg_pos, binsize)
	                res.append("%.4f %.4f %.4f %.4f\n" % (reevx, reevy, reevz, sid))
	                cid[sid] += 1
	#for i in range(bins):
	        #if cid[i] ==0:
	                #cid[i] = 1
	return(cid, countSeg, res)


def main4pp_data(seg, xml, mols, binsize = 0.1, fn=''):
	pos = xml.nodes['position']
	box = xml.box
	bins = int(0.5 * sqrt(box.dot(box))/binsize) + 1
	rg2s = numpy.zeros((bins,), dtype=numpy.float)
	rgn2s = numpy.zeros((bins,),dtype=numpy.float)
	angs = numpy.zeros((bins,),dtype=numpy.float)
	cid = numpy.zeros((bins,),dtype=numpy.float)
	types = xml.nodes['type']
	posA = pos[types=='A']
	cmA = cm_cc(posA, xml.box)
	r_pos = pbc2d(pos-cmA, box)
	#mols = hoomd_mols(xml)
	countSeg = 0
	res = []
	shells = numpy.zeros((bins,))
	for mol in mols.mol_idxes[1:]:
	        length = len(mol)
	        for s in range(0, length, seg):
	                if s + seg > length:
	                        continue
	                countSeg += 1
	                seg_idx = mol[s:s+seg]
	                seg_pos = r_pos[seg_idx]
	                sid = shell_cm(seg_pos, binsize, box)
	                reevx = pbc(seg_pos[-1][0] - seg_pos[0][0], box[0])
	                reevy = pbc(seg_pos[-1][1] - seg_pos[0][1], box[1])
	                reevz = pbc(seg_pos[-1][2] - seg_pos[0][2], box[2])
	                #ss, sdist = shell_dist(seg_pos, binsize, bins)
	                #shells += s
	                #sid = shell_id(seg_pos, binsize)
	                res.append([reevx, reevy, reevz, sid])
	                cid[sid] += 1
	#for i in range(bins):
	        #if cid[i] ==0:
	                #cid[i] = 1
	return(cid, countSeg, res)
