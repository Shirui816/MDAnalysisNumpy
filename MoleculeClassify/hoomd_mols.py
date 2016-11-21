import numpy as np
from Parser.hoomd_xml_pd import hoomd_xml
from MoleculeClassify.bond_hash import bond_hash_unidirect, bond_hash_dualdirect
from MoleculeClassify.body_hash import body_hash_unidirect
from MoleculeClassify.isomer import classify_isomers
from MoleculeClassify.molcules import grabmolecules_without_body, grabmolecules_with_body, grab_iter_dual
from cfunctions.cfunctions.functions import cm_cc

class hoomd_mols(object):
    def bcs(self):
        res = []
        idx = np.arange(self.na)
        for i in self.sbody:
            idxes = idx[self.xml.nodes['body'] == i]
            pos = self.xml.nodes['position'][idxes]
            res.append(cm_cc(pos, self.box))
        return np.array(res)
    
    def __init__(self, hoomd_xml, with_body = False):
        self.xml = hoomd_xml
        self.box = hoomd_xml.box
        self.na = int(hoomd_xml.configure['natoms'].reshape((1,))) # for new feature, the convert from array([1])->1 is deprecated
        self.bond_hash_nn, self.bond_hash_wn = bond_hash_dualdirect(hoomd_xml.nodes['bond'], self.na)
        body = hoomd_xml.nodes.get('body')
        if body is None:
            body= []
        self.body_hash = body_hash_unidirect(body)
        self.sbody = sorted(list(set(body)))
        #gb = grabmolecules_without_body
        #if body!=[]:
            #self.sbody.remove(-1)
            #gb = grabmolecules_with_body
        molecular_list = []
        types = hoomd_xml.nodes.get('type')
        molecular_hash = {}  # save a lot of time ^_^
        idxes = np.arange(self.na)
        bdh = None if not with_body else self.body_hash
        print("Catching Molecules...")
        for i in range(self.na):
            molecular_hash[i] = False
        for i in range(self.na):
            #if molecular_hash[i] == True:
                #continue
            #gb(i, self.body_hash, self.bond_hash_nn, mol_idxes=molecular_idxs,mol_used=molecular_hash)
            molecular_idxs = grab_iter_dual(i, self.bond_hash_nn, mol_used=molecular_hash, body_hash=bdh)#self.body_hash)
            if not len(molecular_idxs)<=1: # avoid monomer and void list
                molecular_list.append(molecular_idxs)
            for atom in molecular_idxs:
                molecular_hash[atom] = True
        #while [] in molecular_list: # this remove operation is really SLOW
            #molecular_list.remove([])
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
