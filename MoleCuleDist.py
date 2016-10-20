import numpy as np
import sys
sys.setrecursionlimit(9500)
from MoleculeClassify.hoomd_mols import hoomd_mols
from Parser.hoomd_xml_pd import hoomd_xml
from sys import argv


class Graph:
	def __init__(self):
		self.neighbors = {}
 
	def add_vertex(self, v):
		if v not in self.neighbors:
			self.neighbors[v] = []
 
	def add_edge(self, u, v):
		self.neighbors[u].append(v)
		# if u == v, do not connect u to itself twice
		#if u != v: # Unconmment if unidirect bond hash was built.
			#self.neighbors[v].append(u)
 
	def vertices(self):
		return list(self.neighbors.keys())
 
	def vertex_neighbors(self, v):
		return self.neighbors[v]

	@staticmethod
	def is_cycl(G):
		Q = []
		V = G.vertices()
		# initially all vertices are unexplored
		layer = { v: -1 for v in V }
		for v in V:
			# v has already been explored; move on
			if layer[v] != -1:
				continue
			# take v as a starting vertex
			layer[v] = 0
			Q.append(v)
			# as long as Q is not empty
			while len(Q) > 0:
				# get the next vertex u of Q that must be looked at
				u = Q.pop(0)
				C = G.vertex_neighbors(u)
				for z in C:
					# if z is being found for the first time
					if layer[z] == -1:
						layer[z] = layer[u] + 1
						Q.append(z)
					elif layer[z] >= layer[u]:
						return True
		return False

def ggm(bond_hash, molecule): # Dual bond_hash
	mol_graph = Graph()
	for atom in molecule:
		mol_graph.add_vertex(atom)
	for atom in molecule:
		for btom in bond_hash[atom]:
			mol_graph.add_edge(atom, btom)
	return mol_graph


from sys import argv

xml = hoomd_xml(argv[1])
mol = hoomd_mols(xml)

def loop_linear(mol):
	loop = []
	linear = []
	for m in mol.mol_idxes:
		m_graph = ggm(mol.bond_hash_nn, m) # Turn molecules into graphs
		if Graph.is_cycl(m_graph): # check if loop in molecule
			loop.append(m)
			#print(m)
		else:
			linear.append(m)

	lloop = [ len(x) for x in loop ]
	llinear = [ len(x) for x in linear ]
	return lloop,llinear
lloop, llinear = loop_linear(mol)

for f in argv[2:]:
	xml = hoomd_xml(f, needed = ['bond','type'])
	mol = hoomd_mols(xml)
	o1, l1 = loop_linear(mol)
	lloop += o1
	llinear += l1


#
## remove the unreacted part
#
#while 5 in llinear:
#	llinear.remove(5)

from pylab import *
#import seaborn as sns
#sns.set(color_codes=True)
f = figure(figsize=(18,9))
binsize = 1
bins_loop = int((max(lloop) - min(lloop))/binsize)
bins_linear = int((max(llinear)-min(llinear))/binsize)
ax1 = f.add_subplot(121)
ax1.hist(lloop, bins=bins_loop, label='Loop')
#sns.distplot(lloop, ax=ax1, label='Loop', rug=1)
ax1.legend()
ax2 = f.add_subplot(122)
ax2.hist(llinear, bins=bins_linear, label='Linear')
#sns.distplot(llinear, ax=ax2, label='Linear', rug=1)
ax2.legend()
ax1.set_xscale('log')
ax2.set_xscale('log')
show()
print(llinear)
