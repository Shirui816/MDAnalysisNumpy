import numpy as np
Bond = {'names': ('name', 'id1', 'id2'), 'formats': ('S16', 'i8', 'i8')}
Angle = {'names': ('name', 'id1', 'id2', 'id3'), 'formats': ('S16', 'i8', 'i8', 'i8')}
Dihedral = {'names': ('name', 'id1', 'id2', 'id3', 'id4'), 'formats': ('S16', 'i8', 'i8', 'i8', 'i8')}
Type = {'names': 'type', 'formats': 'S1'}
Type = np.dtype('S1')
Body = np.dtype('i8')
Mass = np.dtype('f8')
Pos = {'names': ('x', 'y', 'z'), 'formats': ('f8', 'f8', 'f8')}


dtypeDict = dict(bond=Bond, angle=Angle, dihedral=Dihedral, type=Type, body=Body)