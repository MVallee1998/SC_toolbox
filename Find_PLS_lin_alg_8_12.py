import linear_alg_method as lam
import timeit
import numpy as np
import SimplicialComplex as sc
from multiprocessing import Pool

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]


m = 8
n = 4
raw_results_PATH = 'raw_results/PLS_%d_%d' % (m, n)
partial_results_PATH = 'partial_results/PLS_%d_%d' % (m, n)

def text_list(results,path):
    t = open(path, mode='a', encoding='utf-8')
    for K in results:
        t.write(str(K) + '\n')
    t.close()

def text(K,path):
    t = open(path, mode='a', encoding='utf-8')
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


def increment_index_list_type2(index_array,size_index_array):
    k = 0
    while k < size_index_array and index_array[k] == 1:
        index_array[k] = 0
        k += 1
    if k == size_index_array-4:
        print("huitème de passé")
    if k != size_index_array:
        index_array[k] = 1

def new_vect_to_mult_array(vector,size_index_array):
    k=0
    vect_to_mult_array = np.zeros((size_index_array,16384))
    while k<16384 and vector.any()==1:
        vect_to_mult_array[:,k] = vector.copy()
        k+=1
        increment_index_list_type2(vector,size_index_array)
    return vect_to_mult_array

# def new_vect_to_mult_array_1(x,size_index_array,list_2_pow):
#     vect_to_mult_array = np.zeros((size_index_array,4096),dtype=np.float)
#     for k in range(4096):
#         vect_to_mult_array[:, k] = int_to_filter(x,size_index_array,list_2_pow)
#         x+=1
#         if x==2*list_2_pow[size_index_array-1]:
#             x=0
#             break
#     return vect_to_mult_array,x
#
# def int_to_filter(x, nbr_results,list_2_pow):
#     vector = np.zeros(nbr_results,np.int)
#     for k in range(nbr_results):
#         if list_2_pow[k]|x == x:
#             vector[k] = 1
#         if list_2_pow[k]>x:
#             break
#     return vector


if __name__ == '__main__':
    list_char_funct = sc.enumerate_char_funct_orbits(n, m)
    results = []
    eq_classes = []
    for char_funct in list_char_funct:
        start = timeit.default_timer()
        facets = sc.find_facets_compatible_with_lambda(char_funct, m, n)
        M = lam.construct_matrix(facets)
        list_v = lam.find_kernel(M)
        nbr_results = list_v.shape[0]
        print(nbr_results)
        np_facets = np.array(facets)
        A = np.transpose(list_v)
        vect_to_mult = np.zeros(nbr_results, dtype=np.float)
        vect_to_mult[0] = 1
        while vect_to_mult.any() == 1:
            vect_to_mult_array = new_vect_to_mult_array(vect_to_mult, nbr_results)
            candidate_array = A.dot(vect_to_mult_array) % 2
            prod = M.dot(candidate_array)
            having_first_facet = candidate_array[0, :] == 1
            verifying_G_theorem = np.sum(candidate_array, axis=0) <= G_vector[n - 1]
            having_every_closed_ridges = np.logical_not((prod >= 4).any(axis=0))
            good_candidates = candidate_array.T[
                np.logical_and(np.logical_and(having_first_facet, verifying_G_theorem), having_every_closed_ridges)]
            for good_candidate in good_candidates:
                good_candidate_facets = np_facets[good_candidate == 1]
                good_candidate_facets_list = list(good_candidate_facets)
                already_found = False
                K = sc.PureSimplicialComplex(good_candidate_facets_list)
                for data in eq_classes:
                    L, is_PLS, f = data
                    if sc.are_isom(K, L):
                        already_found = True
                        break
                if not already_found:
                    if K.Pic == 4 and K.is_a_seed() and K.is_Z2_homology_sphere() and K.is_promising():
                        eq_classes.append((K, True, char_funct))
                        text(good_candidate_facets_list, partial_results_PATH)
                        results.append(good_candidate_facets_list)
                        print(good_candidate_facets_list, len(results))
                    else:
                        eq_classes.append((K, False, char_funct))
                # K = sc.PureSimplicialComplex(good_candidate_facets_list)
                # text(good_candidate_facets_list,raw_results_PATH)
        stop = timeit.default_timer()
        print("Time spent: ", stop - start)
        print("number of results", len(results))
    text_list(results,raw_results_PATH)

# list_char_funct = sc.enumerate_char_funct_orbits(n, m)
# for char_funct in list_char_funct:
#     results = f(char_funct)
