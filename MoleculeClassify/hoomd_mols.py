import numpy as np
from Parser.hoomd_xml import hoomd_xml
from MoleculeClassify.bond_hash import bond_hash_unidirect
from MoleculeClassify.body_hash import body_hash_unidirect
from MoleculeClassify.isomer import classify_isomers
from MoleculeClassify.molcules import grabmolecules_without_body, grabmolecules_with_body


class hoomd_mols(object):
    def __init__(self, hoomd_xml):
        self.cbox = hoomd_xml.cbox
        self.na = hoomd_xml.configure['natoms'].reshape((1,))
        self.bond_hash_nn, self.bond_hash_wn = bond_hash_unidirect(hoomd_xml.nodes['bond'], self.na)
        self.body_hash = body_hash_unidirect(hoomd_xml.nodes['body'])
        molecular_list = []
        types = hoomd_xml.nodes['type']
        molecular_hash = {}  # save a lot of time ^_^
        idxes = np.arange(self.na)
        print("Catching Molecules...")
        for i in range(self.na):
            molecular_hash[i] = False
        for i in range(self.na):
            if molecular_hash[i] == True:
                continue
            molecular_idxs = []
            #grabmolecules_with_body(i, self.body_hash, self.bond_hash_nn, mol_idxes=molecular_idxs,
                                    #mol_used=molecular_hash)
            grabmolecules_without_body(i, self.bond_hash_nn, mol_idxes=molecular_idxs,mol_used = molecular_hash)
            molecular_list.append(molecular_idxs)
            for atom in molecular_idxs:
                molecular_hash[atom] = True
        while [] in molecular_list:
            molecular_list.remove([])
        self.mol_idxes = np.array([ np.array(x) for x in molecular_list ])
        #print(self.mol_idxes)
        molecular_types = []
        #molecular_bodies = []
        for m in self.mol_idxes:
            molecular_types.append(types[np.array(m)])
            #molecular_bodies.append(self.nodes['body'][np.array(m)])
        print("Done.")
        self.mol_types = np.array([ np.array(x) for x in molecular_types ])
        self.isomer_hash = {}
        print("Classifying isomers...")
        classify_isomers(self.mol_types, self.isomer_hash)
        print("Done.")
