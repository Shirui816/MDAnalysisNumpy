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
cid = loadtxt('cid.txt').T[1]
shells = {}

for i in range(len(cid)):
    if cid[i] == 0:
        cid[i] = 1


print("Loading data...")
data = loadtxt(argv[1])
data = np.append(data, np.zeros(data.shape, dtype=np.complex), axis=0)
print('Data loaded, start computing...')
leng, wid = data.shape
dnframes = int(leng/nvecs)
snframes = int(dnframes/2)

for i in range(cid.shape[0]):
    shells[i] = np.zeros((snframes,), dtype=np.complex)

shell_id = {}
#acf = np.zeros((dnframes,) dtype=np.complex)
acf_dict = {}
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
    acf_dict[i] = a
    shell_id[i] = vec0[0][-1]

for i in range(nvecs):
    sid = shell_id[i]
    acf = acf_dict[i]
    shells[sid] += acf[:snframes]

for i in shells:
    acf = shells[i]/cid[i]
    acf = np.abs(acf)
    acf /= np.arange(snframes, 0, -1)
    acf /= acf[0]
    shells[i] = acf


for i in sorted(shells):
    if i>10:
        break
    o = open('%s.dat' % (i), 'w')
    for t, k in enumerate(shells[i]):
        o.write('%.4f %.4f\n' % (t, k))

print("*** Data saved in ``acf.log'' ***")
print("Done.")
