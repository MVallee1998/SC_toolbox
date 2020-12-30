import linear_alg_method as lam
from itertools import combinations, permutations
import timeit
import numpy as np
import SimplicialComplex as sc
from scipy.sparse import csr_matrix
from numba import jit


# char_funct = [3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15, 8, 4, 2, 1]
# m = 15
# n = 11
# M, facets, ridges = construct_matrix(char_funct, n, m)
# M_sparse = csr_matrix(M)
# print(len(facets), len(ridges))

def mod2(x):
    return x % 2


m = 9
n = 5

def text(result):
    name = 'result/PLS_%d_%d_lin_alg' % (m,n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()
# char_funct = [3, 5, 6, 7, 4, 2, 1]
# m = 7
# n = 4
counter = 0
list_char_funct = sc.enumerate_char_funct_orbits(n,m)
results = []
start = timeit.default_timer()

for char_funct in list_char_funct:
    print(counter/len(list_char_funct)*100)
    counter+=1
    # print(char_funct)
    M, facets, ridges = lam.construct_matrix(char_funct, n, m)
    list_v = lam.find_kernel(M)
    # print(list_v.shape)
    nbr = list_v.shape[0]
    number_results, nbr_facets = list_v.shape
    A = np.transpose(list_v)
    zero = np.zeros(number_results)
    for k in range(1, number_results):
        list_combinations = list(combinations(range(number_results), k))
        # candidates = np.zeros((nbr_facets,len(list_combinations)))
        for l in range(len(list_combinations)):
            linear_comb = list_combinations[l]
            vect_to_mult = np.zeros(number_results)
            vect_to_mult[list(linear_comb)] = np.ones(k)
            candidate = (A.dot(vect_to_mult) % 2)
            prod = M.dot(candidate)
        # if not ((prod >= 4).any(axis = 0)):
        #     if candidate[0] == 1:
        #         info_facets = list(candidate.reshape(nbr_facets))
        #         K = [facets[index] for index in range(len(info_facets)) if info_facets[index] == 1]
        #         if K not in results:
        #             results.append(K)
stop = timeit.default_timer()
print("Time spent: ", stop - start)
print("number of results",len(results))
# text(results)
