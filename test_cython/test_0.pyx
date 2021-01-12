import timeit
import numpy as np
cimport numpy as np
cimport cython


G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_f

@cython.boundscheck(False)  # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function

def increment_index_list_type2(np.ndarray[DTYPE_f, ndim=1] index_array,int size_index_array):
    cdef int k
    k = 0
    while k < size_index_array and index_array[k] == 1:
        index_array[k] = 0
        k += 1
    if k == size_index_array-4:
        print("huitème de passé")
    if k != size_index_array:
        index_array[k] = 1

def new_vect_to_mult_array(np.ndarray[DTYPE_f, ndim=1] vector, int size_index_array):
    cdef int l
    l=0
    cdef np.ndarray vect_to_mult_array
    vect_to_mult_array = np.zeros((size_index_array,2048),np.float)
    while l<2048 and vector.any()==1:
        vect_to_mult_array[:,l] = vector.copy()
        l+=1
        increment_index_list_type2(vector,size_index_array)
    return vect_to_mult_array

def f(np.ndarray[DTYPE_f, ndim=2] M,np.ndarray[DTYPE_f, ndim=1]  np_facets,np.ndarray[DTYPE_f, ndim=2]  list_v,int n):
    cdef double start, stop
    cdef int nbr, r
    cdef np.ndarray[DTYPE_f, ndim=2] A, vect_to_mult_array, candidate_array
    cdef np.ndarray[DTYPE_f, ndim=1] vect_to_mult
    cdef int nbr_candidates
    start = timeit.default_timer()
    number_results = list_v.shape[0]
    print(number_results)
    A = np.transpose(list_v)
    vect_to_mult = np.zeros(number_results,np.float)
    results = []
    vect_to_mult[0] = 1
    while vect_to_mult.any()==1:
        vect_to_mult_array = new_vect_to_mult_array(vect_to_mult, number_results)
        candidate_array = A.dot(vect_to_mult_array) % 2
        prod = M.dot(candidate_array)
        having_first_facet = candidate_array[0,:] == 1
        verifying_G_theorem = np.sum(candidate_array,axis=0) <= G_vector[n - 1]
        having_every_closed_ridges = np.logical_not((prod >= 4).any(axis=0))
        good_candidates = candidate_array.T[np.logical_and(np.logical_and(having_first_facet,verifying_G_theorem),having_every_closed_ridges)]
        nbr_candidates = good_candidates.shape[0]
        for r in range(nbr_candidates):
            good_candidate = good_candidates[r]
            good_candidate_facets = np_facets[good_candidate==1]
            good_candidate_facets_list = list(good_candidate_facets)
            results.append(good_candidate_facets_list)
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    print("number of results", len(results))
    return results