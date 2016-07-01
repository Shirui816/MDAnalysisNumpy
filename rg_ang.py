import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml import hoomd_xml

from sys import argv

fs = argv[1:]
from RgAng import main
xml1 = hoomd_xml(fs[0])
mols = hoomd_mols(xml1)

rg2s, rgn2s, angs, r = main(10, xml1, mols, binsize=0.2)
for fn in fs[1:]:
	xml = hoomd_xml(fn)
	print(fn)
	rg2s1, rgn2s1, angs1, r1 = main(10, xml, mols, binsize=0.2, fn=fn)
	rg2s += rg2s1
	rgn2s += rgn2s1
	angs += angs1
	r += r1





nf = len(argv[1:])
rg2s/=nf
rgn2s/=nf
angs/=nf
r/=nf

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

