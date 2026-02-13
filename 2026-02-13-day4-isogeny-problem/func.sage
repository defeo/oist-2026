def first_nodes(j):
    F = j.parent()
    R.<X> = PolynomialRing(F)
    Phi = classical_modular_polynomial(2)(j, X)
    return Phi.roots(multiplicities=False)

def next_nodes(j0, j1):
    F = j0.parent()
    R.<X> = PolynomialRing(F)
    Phi = classical_modular_polynomial(2)(j1, X)
    assert Phi(j0) == 0
    F = R(Phi/(X - j0))
    return F.roots(multiplicities=False)

def dfs_2isog_all(j0, depth):
    nodes_stack = []
    chain = []
    d = 1
    for j in first_nodes(j0):
        nodes_stack.append((j0, j, d))
    
    ret_chains = []
    while nodes_stack:
        j0, j, dd = nodes_stack.pop()
        for _ in range(d - dd):
            chain.pop()
        d = dd
        if d == depth:
            chain.append(j)
            ret_chains.append(chain[:])
            chain.pop()
            continue
        for jd in next_nodes(j0, j):
            nodes_stack.append((j, jd, d + 1))
        chain.append(j)
    return ret_chains

def dfs_2isog_cond(j0, depth, cond):
    nodes_stack = []
    chain = []
    d = 1
    for j in first_nodes(j0):
        nodes_stack.append((j0, j, d))
    
    ret_chains = []
    while nodes_stack:
        j0, j, dd = nodes_stack.pop()
        for _ in range(d - dd):
            chain.pop()
        d = dd
        if cond(j):
            chain.append(j)
            return chain[:]
        if d == depth:
            continue
        for jd in next_nodes(j0, j):
            nodes_stack.append((j, jd, d + 1))
        chain.append(j)
    return None

def num_of_codomain_j(Fp2, depth):
    j0 = random_ss_j(Fp2)
    chains = dfs_2isog_all(j0, depth)
    S = set([c[-1] for c in chains])
    return len(S)

def random_ss_j(Fp2):
    while True:
        j = Fp2.random_element()
        if EllipticCurve_from_j(j).is_supersingular():
            return j
        
def random_ss_j_from_j0(j0, depth):
    js = first_nodes(j0)
    j1 = js[randint(0, len(js) - 1)]
    for _ in range(depth - 1):
        js = next_nodes(j0, j1)
        j0, j1 = j1, js[randint(0, len(js) - 1)]
    return j1

def exhaustive_search(j1, j2, depth):
    cond = lambda j: j == j2
    return dfs_2isog_cond(j1, depth, cond)

def mitm_search(j1, j2, depth):
    d1 = depth // 2
    d2 = depth - d1
    chains1 = dfs_2isog_all(j1, d1)
    H = dict()
    for c in chains1:
        H[c[-1]] = c
    cond = lambda j: j in H
    chain2 = dfs_2isog_cond(j2, d2, cond)
    if chain2 is None:
        return None
    chain1 = H[chain2[-1]]
    chain2 = chain2[:-1]
    chain2.reverse()
    chain = chain1 + chain2
    Phi = classical_modular_polynomial(2)
    assert Phi(j1, chain[0]) == 0
    assert Phi(chain[-1], j2) == 0
    for k in range(1, len(chain)):
        assert Phi(chain[k - 1], chain[k]) == 0
    return chain

def p_for_cm(prime_list, size):
    while True:
        p = random_prime(size)
        check = True
        for q in prime_list:
            if legendre_symbol(-p, q) != 1:
                check = False
                break
        if check:
            return p

def cm_action(j, prime_list, exp_list):
    Fp = j.parent()
    Fpx.<x> = PolynomialRing(Fp)

    path = [j]
    degs = []
    for q, e in zip(prime_list, exp_list):
        if e > 0:
            Phi = Fpx(classical_modular_polynomial(q)(path[-1], x))
            rs = Phi.roots(multiplicities=False)
            j_new = rs[0]
            path.append(j_new)
            degs.append(q)
            for _ in range(e - 1):
                Phi = Fpx(classical_modular_polynomial(q)(path[-1], x))
                assert Phi(path[-2]) == 0
                Phi = Fpx(Phi/(x - path[-2]))
                rs = Phi.roots(multiplicities=False)
                j_new = rs[0]
                path.append(j_new)
                degs.append(q)
    assert len(path) == len(degs) + 1
    for i in range(len(degs)):
        assert classical_modular_polynomial(degs[i])(path[i], path[i + 1]) == 0
    return path, degs

def GD_search(j1, j2, depth, prime_list, exp_bound):
    Fp2 = j1.parent()
    Fp = Fp2.base_ring()
    p = Fp2.characteristic()
    is_fp = lambda j: j^p == j

    path1 = [j1] + dfs_2isog_cond(j1, depth, is_fp)
    degs1 = [2] * (len(path1) - 1)
    path2 = [j2] + dfs_2isog_cond(j2, depth, is_fp)
    degs2 = [2] * (len(path2) - 1)

    j1_fp = Fp(path1[-1])
    j2_fp = Fp(path2[-1])
    E1 = EllipticCurve_from_j(j1_fp)
    E2 = EllipticCurve_from_j(j2_fp)

    Fpx.<x> = PolynomialRing(Fp)
    f1 = Fpx(E1.defining_polynomial()(x=x, y=0, z=1))
    f2 = Fpx(E2.defining_polynomial()(x=x, y=0, z=1))

    if len(f1.roots()) != 3:
        phi = Fpx(classical_modular_polynomial(2)(j1_fp, x))
        rs = phi.roots(multiplicities=False)
        j1_new = rs[0]
        E1 = EllipticCurve_from_j(j1_new)
        f1 = Fpx(E1.defining_polynomial()(x=x, y=0, z=1))
        assert len(f1.roots()) == 3
        path1.append(j1_new)
        degs1.append(2)
        j1_fp = j1_new
    if len(f2.roots()) != 3:
        phi = Fpx(classical_modular_polynomial(2)(j2_fp, x))
        rs = phi.roots(multiplicities=False)
        j2_new = rs[0]
        E2 = EllipticCurve_from_j(j2_new)
        f2 = Fpx(E2.defining_polynomial()(x=x, y=0, z=1))
        assert len(f2.roots()) == 3
        path2.append(j2_new)
        degs2.append(2)
        j2_fp = j2_new

    for q in prime_list:
        for exp_list in cartesian_product_iterator([range(exp_bound + 1)] * len(prime_list)):
            subpath, subdegs = cm_action(j1_fp, prime_list, exp_list)
            if subpath[-1] == j2_fp:
                path1 += subpath[1:]
                degs1 += subdegs
                path1 += path2[::-1][1:]
                degs1 += degs2[::-1]
                return path1, degs1
