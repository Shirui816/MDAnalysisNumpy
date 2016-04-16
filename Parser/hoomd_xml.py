from numpy import array, loadtxt
import numpy as np
from io import StringIO
import xml.etree.cElementTree as ET  #pypy will be a bit slower than python
import warnings
from DataStructure.DtypeDict import dtypeDict
try:
    from iopro import loadtxt
except ImportError:
    warnings.warn("No module iopro, I can't accelerate while files are large.")
class hoomd_xml(object):
    @staticmethod
    def _get_attrib(dd):
        dt = eval('[' + ','.join(["('%s', int)" % key for key in dd.keys()]) + ']')
        values = [tuple(dd.values())]
        return array(values, dtype=dt)
    def __init__(self, filename):
        tree = ET.ElementTree(file=filename)
        root = tree.getroot()
        configuration = root[0]
        self.configure = self._get_attrib(configuration.attrib)
        self.nodes = {}
        for e in configuration:
            if e.tag == 'box':
                self.box = self._get_attrib(e.attrib)
                self.cbox = np.array([self.box['lx'], self.box['ly'], self.box['lz']]).reshape((3,))
                continue
            try:
                #print(e.tag)
                dt = dtypeDict[e.tag]
                self.nodes[e.tag] = loadtxt(StringIO(e.text), dtype=dt)
            except KeyError:
                self.nodes[e.tag] = loadtxt(StringIO(e.text))
if __name__ == "__main__":
    from sys import argv
    file = argv[1]
    xml = hoomd_xml(argv[1])
    print(xml.nodes['bond']['name'])
    print("----------")
    print(xml.box, xml.configure)
    print(xml.nodes)