import linear_alg_method as lam
import timeit
import numpy as np
import SimplicialComplex as sc
import numba as nb
from multiprocessing import Pool

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]


m = 12
n = 8
raw_results_PATH = 'test_results/PLS_%d_%d' % (m, n)

def text(results,path):
    t = open(path, mode='a', encoding='utf-8')
    for K in results:
        t.write(str(K) + '\n')
    t.close()

@nb.njit
def increment_index_list(arr):
    for idx in range(len(arr)):
        if idx == len(arr)-3:
            print("huitieme")
        if arr[idx]==0:
            arr[idx] = 1
            break
        else:
            arr[idx] = 0

@nb.njit
def new_vect_to_mult_array(vector,size_index_array):
    k=0
    vect_to_mult_array = np.zeros((size_index_array,4096*4))
    while k<4096*4 and vector.any()==1:
        vect_to_mult_array[:,k] = vector.copy()
        k+=1
        increment_index_list(vector)
    return vect_to_mult_array


@nb.njit
def get_product(M,A,vect_to_mult_array):
    candidate_array = np.mod(A.dot(vect_to_mult_array), 2)
    prod = M.dot(candidate_array)
    return candidate_array, prod



def f(char_funct):
    start = timeit.default_timer()
    facets = sc.find_facets_compatible_with_lambda(char_funct,m,n)
    M = lam.construct_matrix(facets)
    list_v = lam.find_kernel(M)
    nbr_results = list_v.shape[0]
    print(nbr_results)
    np_facets = np.array(facets)
    A = np.transpose(list_v)
    vect_to_mult = np.zeros(nbr_results,dtype=np.float)
    results = []
    vect_to_mult[0] = 1
    while vect_to_mult.any()==1:
        vect_to_mult_array = new_vect_to_mult_array(vect_to_mult, nbr_results)
        candidate_array, prod = get_product(M,A, vect_to_mult_array)
        having_first_facet = candidate_array[0, :] == 1
        verifying_G_theorem = np.sum(candidate_array, axis=0) <= G_vector[n - 1]
        having_every_closed_ridges = np.logical_not((prod >= 4).any(axis=0))
        good_candidates = candidate_array.T[
            np.logical_and(np.logical_and(having_first_facet, verifying_G_theorem), having_every_closed_ridges)]
        for good_candidate in good_candidates:
            good_candidate_facets = np_facets[good_candidate==1]
            good_candidate_facets_list = list(good_candidate_facets)
            results.append(good_candidate_facets_list)
            # K = sc.PureSimplicialComplex(good_candidate_facets_list)
            # text(good_candidate_facets_list,raw_results_PATH)
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    print("number of results", len(results))
    return results


# index_array = np.zeros(5)
# index_array[0] = 1
# vect_to_mult_array = new_vect_to_mult_array(index_array,5)
# print(vect_to_mult_array)
# for k in range(10):
#     print(vect_to_mult_array)
#     vect_to_mult_array = new_vect_to_mult_array(index_array, 5)


if __name__ == '__main__':
    list_char_funct = sc.enumerate_char_funct_orbits(n, m)
    with Pool(processes=4) as pool:
        big_result = pool.imap(f, list_char_funct)
        for results in big_result:
            text(results,raw_results_PATH)

# list_char_funct = sc.enumerate_char_funct_orbits(n, m)
# for char_funct in list_char_funct[:1]:
#     results = f(char_funct)
