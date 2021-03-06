import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml

from sys import argv

import os


SEG, BINSIZE = 50, 2


fs = argv[1:]
from Functions.SegAcf import main4pp
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)

cid,cs,res = main4pp(SEG, xml1, mols, binsize=BINSIZE)
alres = {}
alres[0] = res

o = open('cid.txt', 'w')

for i in range(cid.shape[0]):
	o.write("%.4f %.4f\n" % (i, cid[i]))

p = open('seg_num.txt','w')
p.write('%.4f\n' % (cs))



import pp



THS = 5

allfiles = fs[1:]
mm = int(len(allfiles)/THS)
rem = len(allfiles) % THS

job_server = pp.Server(THS)


def run(T):
	resd = {}
	s, e, flist, mols, BINSIZE,SEGMENT = T
	for i in range(s,e):
		fn = flist[i]
		xml = Parser.hoomd_xml_pd.hoomd_xml(fn)
		cid1, cs, res = Functions.SegAcf.main4pp(SEGMENT, xml, mols, binsize=BINSIZE)
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

o = open('ree_pp.txt', 'w')

for k in sorted(alres):
	data = alres[k]
	for i in data:
		o.write(i)
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
