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

try:
    from iopro import loadtxt
    ldtxt = 'iopro.loadtxt'
except ImportError:
    from numpy import loadtxt
    ldtxt = 'numpy.loadtxt'

print()
print("*** You are using fft from %s with up to %s thread(s). Data will be loaded by ``%s()''. ***" % (fft_ver, max_cpu, ldtxt))
print()
from sys import argv
import numpy as np

import sys

if len(argv) != 2:
    print("Usage: python autocorr_vec.py <your-data-file>") 
    sys.exit(0)

nvecs = loadtxt('seg_num.txt')
nvecs = int(nvecs)


print("Loading data...")
data = loadtxt(argv[1])
data = np.append(data, np.zeros(data.shape, dtype=np.complex), axis=0)
print('Data loaded, start computing...')
leng, wid = data.shape
dnframes = int(leng/nvecs)
snframes = int(dnframes/2)


acf = np.zeros((dnframes,), dtype=np.complex)
idx0 = np.arange(0, leng, nvecs)
for i in range(nvecs):
    idx = idx0 + i
    vec0 = data[idx]
    x0 = vec0[:,0]
    y0 = vec0[:,1]
    z0 = vec0[:,2]
    fx = fft(x0)
    fy = fft(y0)
    fz = fft(z0)
    sf = fx * np.conj(fx) + fy * np.conj(fy) + fz * np.conj(fz)
    a = ifft(sf)
    acf += a



acf/=nvecs # Average of spectrums
o = open('acf.log','w')
acf = np.absolute(acf) # Amplitude
nacf = acf[:snframes] / np.arange(snframes, 0, -1)
nacf /= nacf[0]
for i,j in enumerate(nacf):
    o.write("%s %s\n" % (i+1, j))
o.close()
print("*** Data saved in ``acf.log'' ***")
print("Done.")
