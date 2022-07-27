from sympy import *
from sympy import groebner

import SimplicialComplex as sc

x, y, z, t = symbols('x y z t')
from sympy.abc import x, y, z, t
import json

list_gens = [x, y, z, t]


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def construct_ideal(K, IDCM_bin):
    n = K.n
    m = K.m
    p = m - n
    K.compute_MNF_set()
    list_monomials = []
    for k in range(n):
        P = Poly(0 * x, gens=list_gens, domain=GF(2))
        for l in range(p):
            if (IDCM_bin[k] >> l) & 1:
                list_deg = [0] * p
                list_deg[l] = 1
                P += polys.monomials.Monomial(list_deg, gens=list_gens).as_expr()
        list_monomials.append(P.copy())
    list_monomials += list_gens[:p]
    L = []
    for MNF_bin in K.MNF_set_bin:
        P = Poly(1 + x - x, gens=list_gens, domain=GF(2))
        for l in range(m):
            if (MNF_bin >> l) & 1:
                P *= list_monomials[l]
        L.append(P.copy())
    GB = groebner(L, gens=list_gens, domain=GF(2), method='buchberger', order='grevlex')
    return GB

for n_seed in range(2,10):
    results_path = 'final_results/CSPLS_%d_%d' % (n_seed, n_seed + 4)
    list_m_n_seeds = [json.loads(facets_bytes) for facets_bytes in read_file(results_path)]
    for K_facets in list_m_n_seeds:
        K = sc.PureSimplicialComplex(K_facets)
        list_IDCM_bin = sc.IDCM_Garrison_Scott(K)
        if len(list_IDCM_bin) == 1:
            continue
        # print('n=',n_seed)
        list_homology = []
        for IDCM_bin in list_IDCM_bin:
            # print(len(construct_ideal(K, IDCM_bin)))
            H = sc.find_Z4_homology(K,IDCM_bin)
            test=False
            for H0 in list_homology:
                if (H0==H).all():
                    test=True
                    break
            if not test:
                list_homology.append(H)
        if len(list_homology)>1:
            print(n_seed,list_homology)


print('finished')