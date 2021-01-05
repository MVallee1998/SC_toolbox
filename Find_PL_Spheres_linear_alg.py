import linear_alg_method as lam
from itertools import combinations, permutations
import timeit
import numpy as np
import SimplicialComplex as sc
from scipy.sparse import csr_matrix
from multiprocessing import Pool

# char_funct = [3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15, 8, 4, 2, 1]
# m = 15
# n = 11
# M, facets, ridges = construct_matrix(char_funct, n, m)
# M_sparse = csr_matrix(M)
# print(len(facets), len(ridges))
G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]


def mod2(x):
    return x % 2


vmod2 = np.vectorize(mod2)
m = 8
n = 4


def text(result):
    name = 'result/PLS_%d_%d_lin_alg' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


def increment_index_list(index_list, max_size):
    if len(index_list) < max_size:
        if index_list[0] != 0:
            return [0] + index_list
        else:
            k = 0
            while k < len(index_list) and index_list[k] == k:
                k += 1
            return [k] + index_list[k:]
    else:
        return []


def increment_index_list_type2(index_array):
    k = 0
    while k < len(index_array) and index_array[k] == 1:
        index_array[k] = 0
        k += 1
    if k != len(index_array):
        index_array[k] = 1

def increment_index_list_type3(index_array):
    k = 0
    list_zeros = np.where(index_array == 0)[0]
    if list_zeros.size>0:
        # print(index_array)
        first_zero = list_zeros[0]
        index_array[:first_zero] = 0
        index_array[first_zero] = 1
        # print(index_array)
    else:
        index_array[:]=0

def f(char_funct):
    start = timeit.default_timer()
    M, facets, ridges = lam.construct_matrix(char_funct, n, m)
    list_v = lam.find_kernel(M)
    nbr = list_v.shape[0]
    print(nbr)
    number_results, nbr_facets = list_v.shape
    A = np.transpose(list_v)
    vect_to_mult = np.zeros(number_results)
    results = []
    vect_to_mult[0] = 1
    while vect_to_mult.any()==1:
        candidate = A.dot(vect_to_mult) % 2
        prod = M.dot(candidate)
        if candidate[0] == 1:
            if np.where(candidate == 1)[0].size <= G_vector[n - 1]:
                if not ((prod >= 4).any()):
                    info_facets = list(candidate.reshape(nbr_facets))
                    K = [facets[index] for index in range(len(info_facets)) if info_facets[index] == 1]
                    if K not in results:
                        results.append(K)
        increment_index_list_type2(vect_to_mult)
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    print("number of results", len(results))
    return results


# index_array = np.zeros(10)
# index_array[0] = 1
# for k in range(10):
#     print(index_array)
#     increment_index_list_type3(index_array)


if __name__ == '__main__':
    list_char_funct = sc.enumerate_char_funct_orbits(n, m)
    with Pool(processes=6) as pool:
        big_result = pool.imap(f, list_char_funct)
        for results in big_result:
            text(results)
