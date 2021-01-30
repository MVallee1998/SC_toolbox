import numpy as np
from itertools import combinations, permutations
import Z2_linear_algebra
import SimplicialComplex as sc
import timeit
from scipy.sparse import csr_matrix

def construct_matrix(facets):
    ridges = []
    for facet in facets:
        for element in sc.list_2_pow:
            if element | facet == facet:
                ridge = element ^ facet
                if not ridge in ridges:
                    ridges.append(ridge)
    ridges.sort()
    # print("Facets and ridges enumerated")
    M = np.zeros((len(ridges), len(facets)),dtype=np.float)
    for j in range(len(facets)):
        for element in sc.list_2_pow:
            if element | facets[j] == facets[j]:
                ridge = element ^ facets[j]
                i = sc.dichotomie(ridges, ridge)
                M[i, j] = 1
    return M



def create_new_cases_array(v): #not used anymore
    i = np.where(v > 0)[0][-1]
    N = v.shape[0]
    if i < N - 1:
        A = np.vstack([v[:i + 1] for k in range(N - int(i) - 1)])
        return np.vstack((np.transpose(A), np.eye(N - int(i) - 1)))


# create_new_cases_array = np.frompyfunc(create_new_cases_array, 1, 2)


def linear_alg_method(M_sparse, M):  #not used anymore
    def linear_alg_method_rec(array_of_v, results):
        start = timeit.default_timer()
        product = M_sparse.dot(array_of_v)
        stop = timeit.default_timer()
        print("Time spent for product: ", stop - start)
        non_candidate_vectors_filter = (product > 2).any(axis=0)
        candidate_vectors_filter = np.logical_not(non_candidate_vectors_filter)
        pseudo_mfds_filter = np.logical_and(candidate_vectors_filter, np.logical_not((product == 1).any(axis=0)))
        pseudo_mfds = array_of_v[:, pseudo_mfds_filter]
        for pseudo_mfd in np.transpose(pseudo_mfds):
            # print(pseudo_mfd)
            results.append(pseudo_mfd.copy())
        candidate_vectors_filter ^= pseudo_mfds_filter
        start = timeit.default_timer()
        candidate_vectors_filter_index = np.where(candidate_vectors_filter)
        # data = []
        # for index in candidate_vectors_filter_index[0]:
        #     v = array_of_v[:,[index]]
        #     for facet_index in np.where(v==1):
        #         position_unclosed_ridges_of_facet = np.logical_and(M[:,facet_index]==1,product[:,index]==1)
        #         if position_unclosed_ridges_of_facet.any():
        #             for position_one in np.where(position_unclosed_ridges_of_facet)[0]:
        #                 new_facet_index=np.where(M[position_one]==1)[0][0]
        #                 new_v = v.copy()
        #                 new_v[new_facet_index][0] = 1
        #                 # print(v,new_v)
        #                 data.append(new_v)

        # i = np.where(v==[1])[0][-1]

        data = [create_new_cases_array(v) for v in np.transpose(array_of_v[:, candidate_vectors_filter]) if v[-1] != 1]
        stop = timeit.default_timer()
        print("Time spent for creating new cases: ", stop - start)
        if data:
            new_array_of_v = np.hstack(data)
            # print(new_array_of_v)
            linear_alg_method_rec(new_array_of_v, results)

    nbr_ridges, nbr_facets = M_sparse.shape
    original_array_of_v = np.zeros((nbr_facets, 1))
    original_array_of_v[0, 0] = 1
    results = []
    linear_alg_method_rec(original_array_of_v, results)
    print(len(results))

def find_kernel(M):
    nbr_ridges, nbr_facets = M.shape
    N = np.transpose(M).copy()
    i = 0
    j = 0
    basis = np.eye(nbr_facets,dtype=np.float)
    while j < nbr_ridges and i < nbr_facets:
        if N[i, j] == 0:
            i_0 = -1
            for i_iter in range(i + 1, nbr_facets):
                if N[i_iter, j] == 1:
                    i_0 = i_iter
                    break
            if i_0 == -1:
                j += 1
                continue
            N[[i_0, i]] = N[[i, i_0]]
            basis[[i_0, i]] = basis[[i, i_0]]
        for i_iter in range(i + 1, nbr_facets):
            if N[i_iter, j] == 1:
                N[i_iter] = (N[i_iter] + N[i]) % 2
                basis[i_iter] = (basis[i_iter] + basis[i]) % 2
        j += 1
        i += 1
    return basis[(N == 0).all(axis=1)]
