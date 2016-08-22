import numpy as np
import numpy
import pandas as pd

from sys import argv

import os
from sys import exit


SEG = [10, 50,100,200]
SEG = [100]
SLICES = 1


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


from scipy.stats import mode

def acf_shell(cidshape, DATA,fft_ver=fft_ver, max_cpu=max_cpu, SEGNUM=0):
	print()
	print("*** You are using fft from %s with up to %s thread(s). ***" % (fft_ver, max_cpu))
	print()
	nvecs = SEGNUM
	shells = {}
	CID = numpy.zeros(cidshape)
	data = DATA
	data = numpy.append(data, numpy.zeros(data.shape, dtype=numpy.complex), axis=0)
	leng, wid = data.shape
	dnframes = int(leng/nvecs)
	snframes = int(dnframes/2)

	for i in range(cidshape):
		shells[i] = numpy.zeros((snframes,), dtype=numpy.complex)

	shell_id = {}
	acf_dict = {}
	idx0 = numpy.arange(0, leng, nvecs)
	MIDDLE = int(snframes/2)
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
		sids = list((vec0[:,-1].real).astype(numpy.int))[:snframes] # or there will be snframes 0s !!!
		#sid = int(mode(sids)[0][0])
		sid = max(map(lambda val: (sids.count(val), val), set(sids)))[1]
		shell_id[i] = sid # 0 for 1st shell, MIDDLE for the median shell, this is for the most frequent value
		CID[sid] += 1
	

	for i in range(CID.shape[0]):
		if CID[i] == 0:
			CID[i] = 1
	for i in range(nvecs):
		sid = shell_id[i]
		acf = acf_dict[i]
		shells[sid] += acf[:snframes]

	for i in shells:
		acf = shells[i]/CID[i]
		acf = numpy.abs(acf)
		acf /= numpy.arange(snframes, 0, -1)
		acf /= acf[0] if acf[0] != 0 else 1 # avoid NaN
		shells[i] = acf
	return(shells)






def acf_avg(DATA,fft_ver=fft_ver, max_cpu=max_cpu, SEGNUM=0):
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


CIDSHAPE = int(open('cid_shape.txt','r').read())

import sys
for seg in SEG:
	print("Processing on segment %s" % (seg))
	fs_seg = pd.read_csv("%s_ree.txt" % (seg), delim_whitespace=True, squeeze=1, header=None).values
	print("File loaded")
	cs_seg = int(np.loadtxt("%s_seg_num.txt" % (seg)))
	frames_seg = int(fs_seg.shape[0]/cs_seg)
	if frames_seg % SLICES:
		print("Error, %s frames cannot be seperated into %s slices" % (frames_seg, SLICES))
		sys.exit()
	SEP_seg = int(frames_seg/SLICES)
	DATA = fs_seg[0:SEP_seg * cs_seg]
	ACF = acf_avg(fs_seg,fft_ver=fft_ver, max_cpu=max_cpu, SEGNUM=cs_seg) # ACF for all frames, but SLICES for shells.
	SHELLS = acf_shell(CIDSHAPE, DATA, SEGNUM=cs_seg)
	SHELL_COUNT = numpy.zeros(CIDSHAPE)
	for i in SHELLS:
		if SHELLS[i][0] != 0:
			SHELL_COUNT[i] += 1
	for k in range(1, SLICES):
		DATA = fs_seg[SEP_seg * k * cs_seg: SEP_seg * (k+1) * cs_seg]
		#ACF += acf_avg(DATA, fft_ver=fft_ver, max_cpu=max_cpu, SEGNUM=cs_seg)
		shells = acf_shell(CIDSHAPE, DATA, SEGNUM=cs_seg)
		for i in shells:
			if shells[i][0] != 0:
				SHELL_COUNT[i] += 1
			SHELLS[i] += shells[i]
	#ACF /= SLICES
	for i in sorted(SHELLS):
		n = SHELL_COUNT[i] if SHELL_COUNT[i] != 0 else 1
		o = open('%s_%s.dat' % (seg, i), 'w')
		for t, k in enumerate(SHELLS[i]):
			o.write('%.4f %.4f\n' % (t, k/n)) # avoid NaN
		o.close()
	o = open("%s_acf_avg.txt" % (seg), 'w')
	for i,j in enumerate(ACF):
		o.write("%s %s\n" % (i+1, j))
	o.close()

