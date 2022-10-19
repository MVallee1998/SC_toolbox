from sympy import *
from sympy import groebner
import json
import SimplicialComplex as sc
import Z2_linear_algebra as Z2la
import numpy as np

x, y, z, t = symbols('x y z t')
from sympy.abc import x, y, z, t

list_gens = [x, y, z, t]

list_poly_gens = [Poly(x, gens=list_gens, domain=GF(2)),Poly(y, gens=list_gens, domain=GF(2)),Poly(z, gens=list_gens, domain=GF(2)),Poly(t, gens=list_gens, domain=GF(2))]
GL4 = []
for i0 in range(1, 16):
    for i1 in range(1, 16):
        for i2 in range(1, 16):
            for i3 in range(1, 16):
                A_bin = Z2la.Z2Array(4, [i0, i1, i2, i3])
                if A_bin.is_invertible():
                    A = np.zeros((4, 4))
                    A[0] = sc.int_to_bin_array(i0, 4)
                    A[1] = sc.int_to_bin_array(i1, 4)
                    A[2] = sc.int_to_bin_array(i2, 4)
                    A[3] = sc.int_to_bin_array(i3, 4)
                    GL4.append(A.copy())

print(len(GL4))
def transform_base(base, G):
    new_base = []
    for i in range(len(base)):
        for j in range((len(base))):
            if G[i, j] == 1:
                if len(new_base) < i + 1:
                    new_base.append(base[j])
                else:
                    new_base[-1] += base[j]
    return new_base


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def construct_ideal(K, IDCM_bin, test):
    n = K.n
    m = K.m
    p = m - n
    K.compute_MNF_set()
    list_monomials = []
    for k in range(n):
        P = Poly(0 * x, gens=list_gens, domain=GF(2))
        for l in range(p):
            if (IDCM_bin[k] >> l) & 1:
                P += test[p-l-1]
        list_monomials.append(P.copy())
    list_monomials += test[:p]
    # print(list_monomials)
    L = []
    for MNF_bin in K.MNF_set_bin:
        P = Poly(1 + x - x, gens=list_gens, domain=GF(2))
        for l in range(m):
            if (MNF_bin >> l) & 1:
                P *= list_monomials[l]
        L.append(P.copy())
    # print(L)
    GB = groebner(L, gens=list_gens, domain=GF(2), method='buchberger', order='grevlex')
    return GB


for n_seed in range(3, 5):
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
            H = sc.find_Z4_homology(K, IDCM_bin)
            test = False
            for H0 in list_homology:
                if (H0 == H).all():
                    test = True
                    break
            if not test:
                list_homology.append(H)
        if len(list_homology) > 1:
            print("List of Homology classes:", list_homology)
        isom_cohom_P = []
        for i in range(len(list_IDCM_bin)):
            print("global: ", (i / len(list_IDCM_bin)) * 100, '%')
            isom = False
            IDCM_2 = list_IDCM_bin[i]
            for IDCM_1 in isom_cohom_P:
                GB_1 = construct_ideal(K, IDCM_1, list_gens)
                for k in range(len(GL4)):
                    if k % 1000 == 0:
                        print((k / len(GL4)) * 100, '%')
                    G = GL4[k]
                    GB_2 = construct_ideal(K, IDCM_2, transform_base(list_gens, G))
                    if GB_1 == GB_2:
                        print(G)
                        isom = True
                        break
                if isom:
                    break
            if not isom:
                isom_cohom_P.append(IDCM_2.copy())
                if len(isom_cohom_P)>=len(list_homology):
                    print("not an example")
                    break
        print("number of cohomology classes:", len(isom_cohom_P))

print('finished')
