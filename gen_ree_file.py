import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml

from sys import argv

import os

if os.path.isfile('ree.txt'):
	os.remove('ree.txt')


SEG, BINSIZE = 50, 2


fs = argv[1:]
from Functions.SegAcf import main
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)

cid,cs = main(SEG, xml1, mols, binsize=BINSIZE)
for fn in fs[1:]:
	print(fn)
	xml = hoomd_xml(fn)
	cid1, cs = main(SEG, xml, mols, binsize=BINSIZE, fn=fn)
	#cid += cid1

o = open('cid.txt', 'w')

for i in range(cid.shape[0]):
	o.write("%.4f %.4f\n" % (i, cid[i]))

p = open('seg_num.txt','w')
p.write('%.4f\n' % (cs))


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

