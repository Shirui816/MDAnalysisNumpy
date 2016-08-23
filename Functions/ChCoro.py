import numpy as np
from cfunctions.cfunctions.functions import cm_cc, pbc1d, dist1d
import math

def sum_chains(xml, mols, check = 5, thresh=5):
	'''
	Args: xml, mols objects
	Return: 2d-array with data of ith center's jth chain
	'''
	pos = xml.nodes['position']
	mol_indx = mols.mol_idxes
	bcs = mols.bcs()
	box = xml.box
	result = {}
	start = len(bcs)
	result = np.zeros((len(bcs), len(mol_indx) - start))
	for i in range(start, len(mol_indx)):
		mol = mol_indx[i]
		mol_len = len(mol)
		for k in range(len(bcs)):
			bc = bcs[k]
			n = 0
			for j in mol:
				pj = pos[j]
				dist = dist1d(pj, bc, box)
				if dist <= check:
					n += 1
				if n>=thresh:
					result[k][i-start] = 1
					break
	return result
