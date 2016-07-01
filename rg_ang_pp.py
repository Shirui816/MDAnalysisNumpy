import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml

from sys import argv

fs = argv[1:]
from Functions.RgAng import main
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)

SEGMENT,BINSIZE=50,1.0

rg2s, rgn2s, angs, r = main(SEGMENT, xml1, mols, binsize=BINSIZE)
shape = rg2s.shape
#for fn in fs[1:]:
#	xml = hoomd_xml(fn)
#	print(fn)
#	rg2s1, rgn2s1, angs1, r1 = main(10, xml, mols, binsize=0.2, fn=fn)
#	rg2s += rg2s1
#	rgn2s += rgn2s1
#	angs += angs1
#	r += r1

import pp



THS = 30

allfiles = fs[1:]
mm = int(len(allfiles)/THS)
rem = len(allfiles) % THS

job_server = pp.Server(THS)


def run(T):
	s, e, flist, mols, shape, BINSIZE,SEGMENT = T
	rg2st = numpy.zeros(shape)
	rgn2st = numpy.zeros(shape)
	angst = numpy.zeros(shape)
	#rt = numpy.zeros(shape)
	for i in range(s,e):
		fn = flist[i]
		xml = Parser.hoomd_xml.hoomd_xml(fn)
		rg2s1, rgn2s1, angs1, r1 = RgAng.main(SEGMENT, xml, mols, binsize=BINSIZE)
		rg2st += rg2s1
		rgn2st += rgn2s1
		angst += angs1
		#rt += r1
	return(rg2st, rgn2st, angst)

inputs = list((i*mm, i * mm + mm if i<THS-1 else i*mm+rem+mm, allfiles, mols, shape, BINSIZE, SEGMENT) for i in range(THS))
jobs = [(input_, job_server.submit(run,(input_,), modules=("numpy","Parser.hoomd_xml","RgAng"))) for input_ in inputs]



for input_, job in jobs:
	print("LIST",input_[0], input_[1])
	rg2s1, rgn2s1, angs1 = job()
	rg2s += rg2s1
	rgn2s += rgn2s1
	angs += angs1
	#r += r1

job_server.print_stats()








nf = len(argv[1:])
rg2s/=nf
rgn2s/=nf
angs/=nf

rgfile = open('rg.txt', 'w')
rgfile.write('#r, rg2, rgn2\n')
angfile = open('ang.txt', 'w')
angfile.write('#r, ang\n')

for a, b, c in zip(r, rg2s, rgn2s):
	rgfile.write("%.4f %.4f %.4f\n" % (a, b, c))
rgfile.close()

for a, b in zip(r, angs):
	angfile.write("%.4f %.4f\n" % (a,b))
angfile.close()



#mol = hoomd_mols(xml)

#print(mol.__dict__.keys())

#print(mol.isomer_hash.keys())
#print(mol.mol_idxes)
#print(type(mol.mol_idxes[0]),mol.mol_idxes[5],mol.mol_idxes[86])
#print(type(mol.mol_types[0]), mol.mol_types[5],mol.mol_idxes[86])

#pos = xml.nodes['position']
#types = xml.nodes['type']
#posA = pos[types==b'A']
#print(len(posA), posA)

#print(cm(posA, xml.box))

