import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml
from math import sqrt
from sys import argv

import os


SEG, BINSIZE = [10,20,50,100,150,200], 1
for seg in SEG:
	if os.path.isfile("%s_ree.txt" % (seg)):
		os.remove("%s_ree.txt" % (seg))


fs = argv[1:]
from Functions.SegAcf import main
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)

box = xml1.box
bins = int(0.5 * sqrt(box.dot(box))/BINSIZE) + 1

cs = main(SEG, xml1, mols, binsize=BINSIZE)
for fn in fs[1:]:
	print(fn)
	xml = hoomd_xml(fn)
	cs = main(SEG, xml, mols, binsize=BINSIZE, fn=fn)
	#cid += cid1

for seg in SEG:
	p=open('%s_seg_num.txt' % (seg), 'w')
	p.write(str(cs[seg]))
	p.close()

p = open('cid_shape.txt','w')
p.write(str(bins))
p.close()

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

