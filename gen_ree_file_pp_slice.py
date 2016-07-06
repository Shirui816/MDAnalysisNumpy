import numpy as np
import numpy
from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml

from sys import argv

import os
from sys import exit


SEG, BINSIZE = 100, 2
SLICES = 2

fs = argv[1:]

if len(fs) % SLICES:
	print("Cannot slice traj (%s files) into %s pieces!" % (len(fs), SLICES))
	exit(0)
	

from Functions.SegAcf import main4pp_data
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)
ss = int(len(fs)/SLICES)

CIDS = []

for i in range(0, SLICES):
	idx = i + 1
	xmlx = hoomd_xml(fs[i*ss])
	cid, cs, resx = main4pp_data(SEG, xmlx, mols, binsize=BINSIZE)
	CIDS.append(cid)

cid,cs,res = main4pp_data(SEG, xml1, mols, binsize=BINSIZE)
alres = {}
alres[0] = res

SEGNUM = cs


import pp



THS = 20

allfiles = fs[1:]
mm = int(len(allfiles)/THS)
rem = len(allfiles) % THS

job_server = pp.Server(THS)


def run(T):
	resd = {}
	s, e, flist, mols, BINSIZE,SEGMENT = T
	for i in range(s,e):
		fn = flist[i]
		xml = Parser.hoomd_xml_pd.hoomd_xml(fn, needed=['position', 'type'])
		cid1, cs, res = Functions.SegAcf.main4pp_data(SEGMENT, xml, mols, binsize=BINSIZE)
		resd[i+1] = res
	return(resd)

inputs = list((i*mm, i * mm + mm if i<THS-1 else i*mm+rem+mm, allfiles, mols, BINSIZE, SEG) for i in range(THS))
jobs = [(input_, job_server.submit(run,(input_,), modules=("numpy","Parser.hoomd_xml_pd","Functions.SegAcf"))) for input_ in inputs]



for input_, job in jobs:
	print("LIST",input_[0], input_[1])
	resd = job()
	if resd:
		alres.update(resd)

job_server.print_stats()




try:
	from accelerate.mkl.fftpack import rfft, irfft, fft, ifft
	fft_ver = 'mkl'
except ImportError:
	from numpy.fft import fft2,ifft2, ifft, fft
	fft_ver = 'numpy'


try:
	import mkl
	max_cpu = mkl.get_max_threads()
	mkl.set_num_threads(max_cpu)
except ImportError:
	max_cpu = 1
	pass


def acf_shell(cid, DATA,fft_ver=fft_ver, max_cpu=max_cpu, SEGNUM=SEGNUM):
	print()
	print("*** You are using fft from %s with up to %s thread(s). ***" % (fft_ver, max_cpu))
	print()
	nvecs = SEGNUM
	shells = {}

	for i in range(len(cid)):
		if cid[i] == 0:
		    	cid[i] = 1

	data = DATA
	data = numpy.append(data, numpy.zeros(data.shape, dtype=numpy.complex), axis=0)
	leng, wid = data.shape
	dnframes = int(leng/nvecs)
	snframes = int(dnframes/2)

	for i in range(cid.shape[0]):
		shells[i] = numpy.zeros((snframes,), dtype=numpy.complex)

	shell_id = {}
	acf_dict = {}
	idx0 = numpy.arange(0, leng, nvecs)
	for i in range(nvecs):
		idx = idx0 + i
		vec0 = data[idx]
		x0 = vec0[:,0]
		y0 = vec0[:,1]
		z0 = vec0[:,2]
		fx = fft(x0)
		fy = fft(y0)
		fz = fft(z0)
		sf = fx * numpy.conj(fx) + fy * numpy.conj(fy) + fz * numpy.conj(fz)
		a = ifft(sf)
		acf_dict[i] = a
		shell_id[i] = numpy.int(vec0[0][-1].real)

	for i in range(nvecs):
		sid = shell_id[i]
		acf = acf_dict[i]
		shells[sid] += acf[:snframes]

	for i in shells:
		acf = shells[i]/cid[i]
		acf = numpy.abs(acf)
		acf /= numpy.arange(snframes, 0, -1)
		acf /= acf[0]
		shells[i] = acf
	return(shells)






def acf_avg(DATA,fft_ver=fft_ver, max_cpu=max_cpu, SEGNUM=SEGNUM):
	print()
	print("*** You are using fft from %s with up to %s thread(s). ***" % (fft_ver, max_cpu))
	print()
	nvecs = SEGNUM

	data = DATA
	data = numpy.append(data, numpy.zeros(data.shape, dtype=numpy.complex), axis=0)
	leng, wid = data.shape
	dnframes = int(leng/nvecs)
	snframes = int(dnframes/2)
	acf = np.zeros((dnframes,), dtype=np.complex)
	idx0 = numpy.arange(0, leng, nvecs)
	for i in range(nvecs):
		idx = idx0 + i
		vec0 = data[idx]
		x0 = vec0[:,0]
		y0 = vec0[:,1]
		z0 = vec0[:,2]
		fx = fft(x0)
		fy = fft(y0)
		fz = fft(z0)
		sf = fx * numpy.conj(fx) + fy * numpy.conj(fy) + fz * numpy.conj(fz)
		a = ifft(sf)
		acf += a
	acf/=nvecs # Average of spectrums
	o = open('acf.log','w')
	acf = numpy.absolute(acf) # Amplitude
	nacf = acf[:snframes] / numpy.arange(snframes, 0, -1)
	nacf /= nacf[0]
	return(nacf)




SEP = int(len(fs) / SLICES)
L = sorted(alres)
idx = 0
cid = CIDS[0]
DATA = []
for i in range(0, 0+SEP):
        DATA.extend(alres[L[i]])
DATA = numpy.array(DATA, dtype=numpy.complex)
SHELLS = acf_shell(cid, DATA)
acf = acf_avg(DATA)


for k in range(SEP, len(fs), SEP):
	idx = int(k/SEP)
	cid = CIDS[idx]
	DATA = []
	for i in range(k, k+SEP):
		DATA.extend(alres[L[i]])
	DATA = numpy.array(DATA, dtype=numpy.complex)
	shells = acf_shell(cid, DATA)
	for i in shells:
		SHELLS[i] += shells[i]
	acf += acf_avg(DATA)


for i in sorted(SHELLS):
	if i>10:
		break
	o = open('%s.dat' % (i), 'w')
	for t, k in enumerate(SHELLS[i]):
		o.write('%.4f %.4f\n' % (t, k/SLICES))
	o.close()
o = open('acf.log','w')
for i,j in enumerate(acf):
	o.write("%s %s\n" % (i+1, j/SLICES))
o.close()

#mol = hoomd_mols(xml)

#print(mol.__dict__.keys())

#print(mol.isomer_hash.keys())
#print(mol.mol_idxes)
#print(type(mol.mol_idxes[0]),mol.mol_idxes[5],mol.mol_idxes[86])
#print(type(mol.mol_types[0]), mol.mol_types[5],mol.mol_idxes[86])

#from Functions.cm_with_pbc import cm
#pos = xml.nodes['position']
#types = xml.nodes['type']
#posA = pos[types==b'A']
#print(len(posA), posA)

#print(cm(posA, xml.box))
