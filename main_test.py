import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml import hoomd_xml

from sys import argv

xml = hoomd_xml(argv[1])
mol = hoomd_mols(xml)

print(mol.__dict__.keys())

print(mol.isomer_hash.keys())
print(mol.mol_idxes)
print(type(mol.mol_idxes[0]),mol.mol_idxes[5],mol.mol_idxes[86])
print(type(mol.mol_types[0]), mol.mol_types[5],mol.mol_idxes[86])

from Functions.cm_with_pbc import cm
pos = xml.nodes['position']
types = xml.nodes['type']
posA = pos[types==b'A']
print(len(posA), posA)

print(cm(posA, xml.cbox))