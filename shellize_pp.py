from Functions.shellize import shells_c0
from Functions.pbc import pbc1d
from cfunctions.functions import pbc2d
import numpy as np
from Parser.hoomd_xml import hoomd_xml
from sys import argv
import pp
from cfunctions.functions import cm_cc


shells = 200

xml = hoomd_xml(argv[1])
pos = xml.nodes['position']
types = xml.nodes['type']
posC = pos[types == 'C']
m, n = posC.shape
result = np.zeros((m, shells), dtype=np.int)


THS = 5

allfiles = argv[1:]
mm = int(len(allfiles)/THS)
rem = len(allfiles) % THS

job_server = pp.Server(THS)


def run(T):
    shells, m, s, e, flist = T
    res = numpy.zeros((m, shells), dtype=numpy.int)
    for i in range(s,e):
        f = flist[i]
        xml = Parser.hoomd_xml.hoomd_xml(f)
        pos = xml.nodes['position']
        types = xml.nodes['type']
        posA = pos[types=='A']
        center = cm_cc(posA, xml.box)
        posC = pos[types == 'C']
        posC_c0 = pbc2d(posC - center, xml.box)
        res += shells_c0(posC_c0, xml.box, s = shells)
    return(res)


inputs = list((shells, m, i*mm, i * mm + mm if i<THS-1 else i*mm+rem+mm, allfiles) for i in range(THS))
jobs = [(input_, job_server.submit(run,(input_,),depfuncs=(cm_cc, shells_c0, pbc2d, pbc1d), modules=("numpy","Parser.hoomd_xml"))) for input_ in inputs]



for input_, job in jobs:
    print("LIST",input_[2], input_[3], argv[1:][input_[2]:input_[3]])
    result += job()
job_server.print_stats()

np.savetxt('shells.dat',result,fmt='%d')

