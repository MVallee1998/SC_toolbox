import SimplicialComplex as sc
import linear_alg_method as lam
import numpy as np
import test_0

m=9
n=5

def text(result):
    name = 'tests/PLS_%d_%d_lin_alg' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()

list_char_funct = sc.enumerate_char_funct_orbits(n, m)
for char_funct in list_char_funct:
    M, facets, ridges = lam.construct_matrix(char_funct, n, m)
    np_facets = np.array(facets,dtype = np.float)
    list_v = lam.find_kernel(M)
    results = test_0.f(M, np_facets, list_v, n)