from Functions.shellize import shells_c0
from Functions.cm_with_pbc import cm
from Functions.pbc import pbc2d, pbc1d
import numpy as np
from Parser.hoomd_xml import hoomd_xml
from sys import argv


shells = 200

xml = hoomd_xml(argv[1])
pos = xml.nodes['position']
types = xml.nodes['type']
posC = pos[types == b'C']
m, n = posC.shape
result = np.zeros((m, shells), dtype=np.int)




for f in argv[1:]:
    xml = hoomd_xml(f)
    pos = xml.nodes['position']
    types = xml.nodes['type']
    posA = pos[types==b'A']
    center = cm(posA, xml.box)
    posC = pos[types == b'C']
    posC_c0 = pbc2d(posC - center, xml.box)
    result += shells_c0(posC_c0,xml.box, s=shells)


np.savetxt('shells_nopp.dat', result, fmt='%d')