import numpy as np

from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml import hoomd_xml
from Functions.cm_with_pbc import cm, cm_c
from sys import argv
from cfunctions.cfunctions.functions import cm_cc

xml=  hoomd_xml(argv[1])
mol = hoomd_mols(xml)

print(mol.__dict__.keys())

print(mol.isomer_hash.keys())
print(mol.mol_idxes)
print(type(mol.mol_idxes[0]),mol.mol_idxes[5],mol.mol_idxes[86])
print(type(mol.mol_types[0]), mol.mol_types[5],mol.mol_idxes[86])

pos = xml.nodes['position']
types = xml.nodes['type']
posA = pos[types=='A']
print(len(posA), posA)

print(cm(posA, xml.box))
print(cm_c(posA, xml.box))
print(cm_cc(posA, xml.box))

