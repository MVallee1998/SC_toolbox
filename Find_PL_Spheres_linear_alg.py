import linear_alg_method as lam
import timeit
import numpy as np
import cupy as cp
import SimplicialComplex as sc
import numba as nb
from multiprocessing import Pool

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]

np_arrange = np.arange(0,256)
np_arrange_odd = 2*np.arange(0,127) + 1
m = 11
n = 7
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
    vect_to_mult_array = np.zeros((size_index_array,4096*16))
    while k<4096*16 and vector.any()==1:
        vect_to_mult_array[:,k] = vector.copy()
        k+=1
        increment_index_list(vector)
    return vect_to_mult_array


def increment_vect(vect,size_kernel):
    vect[0]+=1
    k=0
    index_maxi = size_kernel // 8 - 2
    pow_maxi = size_kernel % 8
    while k<=index_maxi and vect[k]==0:
        k+=1
        vect[k]+=1
    if k == index_maxi and vect[k] == sc.list_2_pow[pow_maxi]:
        return True
    return False


def new_vect_to_mult_array_1(vect,size_kernel):
    combinations = np.array(np.meshgrid(np_arrange_odd, np_arrange,vect[0],vect[1],vect[2],vect[3],vect[4],vect[5]),dtype = np.uint8).T.reshape(-1,8)
    result = np.unpackbits(combinations,axis=1,bitorder='little', count= size_kernel).T
    if increment_vect(vect,size_kernel):
        return result, False
    return result, True



# @nb.njit
def get_product(M,A,vect_to_mult_array):
    candidate_array = A.dot(vect_to_mult_array) % 2
    prod = M.dot(candidate_array)
    return candidate_array, prod



def f(char_funct):
    start = timeit.default_timer()
    facets = sc.find_facets_compatible_with_lambda(char_funct,m,n)
    M = lam.construct_matrix(facets)
    list_v = lam.find_kernel(M)
    M_cp = cp.asarray(M)
    nbr_results = list_v.shape[0]
    print(nbr_results)
    np_facets = cp.array(facets)
    A = cp.asarray(np.transpose(list_v))
    vect_to_mult = np.zeros(nbr_results)
    results = []
    vect = np.array(np.zeros((6,1)),dtype=np.uint8)
    # vect_to_mult[0] = 1
    # while vect_to_mult.any()==1:
    keep_going = True
    while keep_going:
        # vect_to_mult_array = new_vect_to_mult_array(vect_to_mult, nbr_results)
        vect_to_mult_array, keep_going = new_vect_to_mult_array_1(vect,nbr_results)
        candidate_array, prod = get_product(M_cp,A, cp.asarray(vect_to_mult_array))
        verifying_G_theorem = cp.sum(candidate_array, axis=0) <= G_vector[n - 1]
        having_every_closed_ridges = cp.logical_not((prod >= 4).any(axis=0))
        good_conditions = cp.flatnonzero(cp.logical_and(verifying_G_theorem, having_every_closed_ridges))
        good_candidates = candidate_array[:,good_conditions].T
        for good_candidate in good_candidates:
            good_candidate_facets = np_facets[good_candidate==1]
            good_candidate_facets_list = good_candidate_facets.tolist()
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


# if __name__ == '__main__':
#     list_char_funct = sc.enumerate_char_funct_orbits(n, m)
#     with Pool(processes=6) as pool:
#         big_result = pool.imap(f, list_char_funct)
#         for results in big_result:
#             text(results,raw_results_PATH)

list_char_funct = sc.enumerate_char_funct_orbits(n, m)
for char_funct in list_char_funct[:1]:
    results = f(char_funct)
    text(results,raw_results_PATH)
