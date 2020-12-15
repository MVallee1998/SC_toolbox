import numpy as np
import cupy as cp
from itertools import combinations, permutations
import Z2_linear_algebra
import SimplicialComplex as sc
import timeit
from scipy.sparse import csr_matrix


def construct_matrix(char_function, n, m):
    Pic = m - n
    cofacets = []
    for cofacet_iter in combinations(range(1, m + 1), Pic):
        sub_array = []
        for index in cofacet_iter:
            sub_array.append(char_function[index - 1])
        if Z2_linear_algebra.Z2Array(Pic, sub_array.copy()).is_invertible():
            cofacets.append(list(cofacet_iter))
    facets = []
    list_2_pow = [1]
    for k in range(1, m + 1):
        list_2_pow.append(list_2_pow[-1] * 2)
    for cofacet in cofacets:
        facets.append((list_2_pow[m] - 1) ^ sc.face_to_binary(cofacet, m))
    facets.sort()
    ridges = []
    for facet in facets:
        for element in list_2_pow:
            if element | facet == facet:
                ridge = element ^ facet
                if not ridge in ridges:
                    ridges.append(ridge)
    ridges.sort()
    print("Facets and ridges enumerated")
    M = np.zeros((len(ridges), len(facets)))
    for j in range(len(facets)):
        for element in list_2_pow:
            if element | facets[j] == facets[j]:
                ridge = element ^ facets[j]
                i = sc.dichotomie(ridges, ridge)
                M[i, j] = 1
    return M, facets, ridges





def create_new_cases_array(v):
    i = np.where(v > 0)[0][-1]
    N = v.shape[0]
    if i < N -1:
        A = np.vstack([v[:i + 1] for k in range(N - int(i) - 1)])
        return np.vstack((np.transpose(A), np.eye(N - int(i) - 1)))


# create_new_cases_array = np.frompyfunc(create_new_cases_array, 1, 2)


def linear_alg_method(M_sparse,M):
    def linear_alg_method_rec(array_of_v, results):
        start = timeit.default_timer()
        product = M_sparse.dot(array_of_v)
        stop = timeit.default_timer()
        print("Time spent for product: ", stop - start)
        non_candidate_vectors_filter = (product > 2).any(axis=0)
        candidate_vectors_filter = np.logical_not(non_candidate_vectors_filter)
        pseudo_mfds_filter = np.logical_and(candidate_vectors_filter,np.logical_not((product == 1).any(axis=0)))
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

        data = [create_new_cases_array(v) for v in np.transpose(array_of_v[:, candidate_vectors_filter]) if v[-1]!=1]
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
    print(nbr_ridges,nbr_facets)
    N = np.transpose(M).copy()
    i=0
    j=0
    basis = np.eye(nbr_facets)
    while j<nbr_ridges:
        if N[i,j]==0:
            i_0 = -1
            for i_iter in range(i+1, nbr_facets):
                if N[i_iter,j]==1:
                    i_0 = i_iter
                    break
            if i_0 == -1:
                j+=1
                continue
            N[[i_0,i]] = N[[i,i_0]]
            basis[[i_0, i]] = basis[[i, i_0]]
        for i_iter in range(i+1,nbr_facets):
            if N[i_iter,j]==1:
                N[i_iter] = (N[i_iter]+ N[i]) % 2
                basis[i_iter] = (basis[i_iter]+ basis[i])%2
        j+=1
        i+=1
    print(np.all(M.dot(np.transpose(basis))%2 == np.transpose(N)))
    return basis[(N==0).all(axis=1)]






char_funct = [3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15, 8, 4, 2, 1]
m = 15
n = 11
# M, facets, ridges = construct_matrix(char_funct, n, m)
# M_sparse = csr_matrix(M)
# print(len(facets), len(ridges))




# char_funct = [3, 5, 6, 7, 4, 2, 1]
# m = 7
# n = 4



# facets, ridges = sc.enumerate_facets_and_ridges(char_funct,n,m)
# for i in range(8,len(facets)):
#     print(i)
#     possibilities = combinations(facets,i)
#     for K_iter in possibilities:
#         K = sc.PureSimplicialComplex(list(K_iter))
#         list_lambdas = sc.Garrison_Scott(K)
#         if len(list_lambdas) == 1:
#             print(K_iter)

M, facets, ridges = construct_matrix(char_funct, n, m)
# M_sparse = csr_matrix(M)
print(len(facets), len(ridges))


# M = np.array([[1,1,0,0],[1,0,1,0],[0,1,1,0],[1,0,0,1],[0,1,0,1],[0,0,1,1]])
list_v = find_kernel(M)
print(list_v.shape)
number_results,nbr_facets = list_v.shape
A = np.transpose(list_v)
zero = np.zeros(number_results)
results = []
for k in range(number_results):
    print(np.where(M.dot(list_v[k])==4)[0].size)


# for k in range(6,number_results):
#     print(k)
#     print(len(list(combinations(range(number_results),k))))
#     start = timeit.default_timer()
#     for linear_comb in combinations(range(number_results),k):
#         vector_mult = zero.copy()
#         vector_mult[list(linear_comb)] = np.ones(k)
#         candidate = A.dot(vector_mult.reshape(number_results,1)) % 2
#         prod = M.dot(candidate)
#         if not ((prod==4).any()):
#             if candidate[0]==1:
#                 # print(len(np.where(candidate==1)[0]))
#                 results.append(prod.copy())
#     stop = timeit.default_timer()
#     print("Time spent: ", stop - start)
# print(len(results))
# print(M_sparse)

# linear_alg_method(M_sparse,M)

# for k in range(1):
#     v = cp.random.randint(2, size=(len(facets), 200000))
#     product = M_sparse.dot(v)
#     candidate_vectors_filter = (product <= 2).all(axis=0)


# v = np.array([[1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]])
# print(create_new_cases_array(v, 0))
