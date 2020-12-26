import linear_alg_method as lam
from itertools import combinations, permutations
import timeit
import numpy as np
import SimplicialComplex as sc

# char_funct = [3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15, 8, 4, 2, 1]
# m = 15
# n = 11
# M, facets, ridges = construct_matrix(char_funct, n, m)
# M_sparse = csr_matrix(M)
# print(len(facets), len(ridges))

m = 10
n = 6

def text(result):
    name = 'result/PLS_%d_%d_lin_alg' % (m,n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()
# char_funct = [3, 5, 6, 7, 4, 2, 1]
# m = 7
# n = 4
mini = 100
maxi = 0
list_char_funct = sc.enumerate_char_funct_orbits(n,m)
results = []
for char_funct in list_char_funct:
    start_sub = timeit.default_timer()
    # for perm_iter in permutations([3, 5, 6, 7]):
    #     tests.append(list(perm_iter) + [4, 2, 1])
    # print(char_funct)
    start = timeit.default_timer()
    print(char_funct)
    M, facets, ridges = lam.construct_matrix(char_funct, n, m)

    list_v = lam.find_kernel(M)
    # print(list_v.shape)
    nbr = list_v.shape[0]
    number_results, nbr_facets = list_v.shape
    A = np.transpose(list_v)
    zero = np.zeros(number_results)

    for k in range(1, number_results):
        print(k)
        list_combinations = list(combinations(range(number_results), k))
        for linear_comb in list_combinations:
            candidate = np.zeros((nbr_facets,1))
            for k in linear_comb:
                candidate = (candidate + A[:,k].reshape(nbr_facets,1)) % 2
            prod = M.dot(candidate)
            if not ((prod >= 4).any()):
                if candidate[0] == 1:
                    info_facets = list(candidate.reshape(nbr_facets))
                    K = [facets[index] for index in range(len(info_facets)) if info_facets[index] == 1]
                    if K not in results:
                        results.append(K)
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    print("number of results",len(results))
text(results)
