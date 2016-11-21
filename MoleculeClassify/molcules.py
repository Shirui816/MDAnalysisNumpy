


def grabmolecules_with_body(i, body_hash, bond_hash, mol_idxes, mol_used):
    if mol_used[i] == True:
        return None
    mol_used[i] = True
    mol_idxes.append(i)
    for j in body_hash[i]:
        mol_idxes.append(j)
        for k in bond_hash[j]:
            grabmolecules_with_body(k, body_hash, bond_hash, mol_idxes, mol_used)
    for j in bond_hash[i]:
        grabmolecules_with_body(j, body_hash, bond_hash, mol_idxes, mol_used)

def grabmolecules_without_body(i, body_hash, bond_hash, mol_idxes,  mol_used):
    #print(type(mol_used), 'typeof mol_used', i)
    #print(mol_used[i])
    if mol_used[i] == True:
        return None
    mol_used[i] = True
    mol_idxes.append(i)
    # mol_bodies.append(bodies[i])
    for j in bond_hash[i]:
        grabmolecules_without_body(j,body_hash, bond_hash, mol_idxes, mol_used)

def grab_iter_dual(i, bond_hash, mol_used, body_hash=None):
    S = [i]
    R = []
    while S:
        v = S.pop()
        if mol_used[v] == False:
            R.append(v)
            mol_used[v] = True
            #for w in bond_hash[v]:
                #S.append(w)
            S.extend(bond_hash[v])
            if not body_hash:
                continue
            for w in body_hash[v]:
                S.append(w)
                for x in bond_hash[w]:
                    S.append(x)
    return R

def grab_iter_uni(i, bond_hash, body_hash):
    seen = {i}
    stack = [i]
    result = []
    if not body_hash:
        body_hash = {}
    while stack:
        x = stack[-1]
        nxt = set(bond_hash.get(x, []))|set(body_hash.get(x, [])) - seen
        if nxt:
            #stack.extend(list(nxt))
            stack.extend(nx for nx in bond_hash[x] if nx in nxt) # Ordered
            seen |= nxt
        else:
            result.append(x)
            stack.pop()
    return result
