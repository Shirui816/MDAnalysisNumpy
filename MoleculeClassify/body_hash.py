import numpy as np

def body_hash_unidirect(body):
    body_hash = {}
    print('Building body hash...')
    if body == []:
        print('Done.')
        return None
    natoms = len(body)
    idxes = np.arange(natoms)
    if body == []:
        bodies = []
    else:
        bodies = list(set(list(body[body!=-1])))
    for i in range(natoms):
        body_hash[i] = []
    for b in bodies:
        ids = list(idxes[body == b])
        body_hash[ids[0]]+=ids[1:]
    print('Done.')
    return body_hash

def body_hash_new(body):
    body_hash = {}
    print('Build body hash...')
    natoms = len(body)
    bodies = list(set(list(body)))
    bodies.remove(-1)
    idxes = np.arange(natoms)
    for b in bodies:
        body_hash[b] = idxes[body==b]
    print('Done.')
    return body_hash
