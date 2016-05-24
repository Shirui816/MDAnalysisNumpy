import numpy as np
from cm_with_pbc.py import cm
from pbc import pbc2d, pbc1d


def shells(center, pos, box, s = 300):
	n, m = pos.shape # Ensures that pos is a 2D array
	result = np.zeros((n, s))
	L = box.min()/2
	delta = L/s
	pbc_pos = pbc2d(pos - center, box)
	for i in range(n):
		d = sum(pbc_pos[i]**2)**0.5
		if not d > L:
			result[i][int(d/delta)] += 1
	return(result)

def shells_c0(pos_c0, box, s = 300):
	n, m = pos_c0.shape # Ensures that pos is a 2D array
	result = np.zeros((n, s))
	L = box.min()/2
	delta = L/s
	for i in range(n):
		d = sum(pos_c0[i]**2)**0.5
		if not d > L:
			result[i][int(d/delta)] += 1
	return(result)