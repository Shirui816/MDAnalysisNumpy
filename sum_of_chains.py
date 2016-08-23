
import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml

from sys import argv

import os




fs = argv[1:]
from Functions.ChCoro import sum_chains
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)
bcs = len(mols.bcs())
allmol = len(mols.mol_idxes)
matrix = allmol - bcs


res = [sum_chains(xml1, mols, check = 5, thresh = 5)]

for f in fs[1:]:
	print(f)
	xml = hoomd_xml(f, needed=['position'])
	r = sum_chains(xml, mols, check=5, thresh =5)
	res.append(r)

ts = len(res)
allacf = np.zeros((ts-1,))
normalize = np.arange(ts-1, 0, -1)
for i in range(matrix):
	for j in range(bcs):
		stats = [ M[j][i] for M in res ]
		acf = np.correlate(stats, stats, 'full')[ts:]
		allacf += acf
allacf /= matrix * bcs * normalize

o = open('contact_acf.txt', 'w')
for i, j in enumerate(allacf):
	o.write("%.4f %.4f\n" % (i,j))
o.close()





