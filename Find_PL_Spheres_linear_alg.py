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


# char_funct = [7, 11, 13, 14, 15, 8, 4, 2, 1]
# m = 9
# n = 5

char_funct = [3, 5, 6, 7, 4, 2, 1]
m = 7
n = 4


start = timeit.default_timer()
M, facets, ridges = lam.construct_matrix(char_funct, n, m)

list_v = lam.find_kernel(M)
print(list_v.shape)
number_results, nbr_facets = list_v.shape
A = np.transpose(list_v)
zero = np.zeros(number_results)
results = []

for k in range(1, number_results):
    print(k)
    list_combinations = list(combinations(range(number_results), k))
    for linear_comb in list_combinations:
        candidate = np.zeros((nbr_facets,1))
        for k in linear_comb:
            candidate = (candidate + A[:,k].reshape(nbr_facets,1)) % 2
        prod = M.dot(candidate)
        if not ((prod == 4).any()):
            if candidate[0] == 1:
                info_facets = list(candidate.reshape(nbr_facets))
                K = sc.PureSimplicialComplex(
                    [facets[index] for index in range(len(info_facets)) if info_facets[index] == 1])
                if K.Pic == 3 and K.is_promising() and K.is_Z2_homology_sphere():
                    results.append(K)

stop = timeit.default_timer()
print("Time spent: ", stop - start)
print(len(results))
