import linear_alg_method as lam
import Puzzle_algorithm as puzalg
import json
import numpy as np
import numpy as cp
import SimplicialComplex as sc
import numba as nb
from itertools import combinations
import sys


pow_2 = np.ones(64,dtype=np.uint)
for k in range(62,-1,-1):
    pow_2[k] = pow_2[k+1]*2
# cfoo(pow_2)
np.set_printoptions(threshold=sys.maxsize)

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]

np_arrange = np.arange(0, 256)
np_arrange_odd = 2 * np.arange(0, 127) + 1
m = 15
n = 11
p=m-n
number_steps = 1

raw_results_PATH = 'raw_results/PLS_%d_%d_non_seeds_base' % (m, n)


def text(results, path):
    t = open(path, mode='a', encoding='utf-8')
    for K in results:
        t.write(str(K) + '\n')
    t.close()

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data



def get_product(M, A, vect_to_mult_array):
    candidate_array = cp.mod(A.dot(vect_to_mult_array), 2)
    prod = M.dot(candidate_array)
    return candidate_array, prod


def Gauss(M):
    N = M.copy()
    a, b = M.shape
    i = 0
    j = 0
    while i < a and j < b:
        if N[i, j] == 0:
            i_0 = -1
            for i_iter in range(i + 1, a):
                if N[i_iter, j] == 1:
                    i_0 = i_iter
                    break
            if i_0 == -1:
                j += 1
                continue
            N[[i_0, i]] = N[[i, i_0]]
        for i_iter in range(0, i):
            N[i_iter] = (N[i_iter] + N[i] * N[i_iter, j]) % 2
        for i_iter in range(i + 1, a):
            if N[i_iter, j] == 1:
                N[i_iter] = (N[i_iter] + N[i]) % 2
        j += 1
        i += 1
    return (N)


def reduce_wrt_columns(M, array_columns, starting_row):
    a, c = M.shape
    b = array_columns.size
    i = starting_row
    j = 0
    while i < a and j < b:
        if M[i, array_columns[j]] == 0:
            i_0 = -1
            for i_iter in range(i + 1, a):
                if M[i_iter, array_columns[j]] == 1:
                    i_0 = i_iter
                    break
            if i_0 == -1:
                j += 1
                continue
            M[[i_0, i]] = M[[i, i_0]]
        for i_iter in range(i):
            M[i_iter] = (M[i_iter] + M[i] * M[i_iter, array_columns[j]]) % 2
        for i_iter in range(i + 1, a):
            M[i_iter] = (M[i_iter] + M[i] * M[i_iter, array_columns[j]]) % 2
        j += 1
        i += 1


def rank(M):
    return (np.count_nonzero(np.sum(Gauss(M), axis=1)))


@nb.njit
def give_next_vect(vect, base):
    index = 0
    vect[index] = (vect[index] + 1) % base[index]
    while vect[index] == 0 and index < vect.size - 1:
        index += 1
        vect[index] = (vect[index] + 1) % base[index]
    if index == vect.size - 1:
        if vect[index] == 0:
            return False
    return True


def partition_generators(list_v_new, starting_row, list_not_together, list_distinct_not_together, sum_of_not_together,
                         list_gener_not_together):
    find_new_one = True
    while find_new_one:
        find_new_one = False
        for not_together in list_not_together:
            index_not_together = np.flatnonzero(not_together)
            if np.sum(np.multiply(not_together, sum_of_not_together)) == 0 and rank(
                    list_v_new[starting_row:, index_not_together]) == index_not_together.size - 1:
                sum_of_not_together += not_together
                reduce_wrt_columns(list_v_new, index_not_together, starting_row)
                list_gener_not_together.append(list(range(starting_row, starting_row + index_not_together.size - 1)))
                starting_row += index_not_together.size - 1
                list_distinct_not_together.append(not_together.copy())
                find_new_one = True
                break
    return starting_row


def new_f(K_max_facets,K_sus_facets,list_results):
    M = lam.construct_matrix(K_max_facets)
    nbr_ridges, nbr_facets = M.shape
    np_K_max_facets = np.array(K_max_facets)
    list_v = lam.find_kernel(M.copy())
    nbr_generators = list_v.shape[0]
    I = []
    J = []
    I_cup_J = []
    for i in range(len(K_max_facets)):
        facet = K_max_facets[i]
        if 1&facet:
            I_cup_J.append(i)
            if facet in K_sus_facets:
                I.append(i)
            else:
                J.append(i)
    reduce_wrt_columns(list_v, np.array(I+J), 0)
    # print(list_v[:,J].T)
    if (np.mod(np.sum(list_v[:,I],axis=0).T,2)!=1).any() or (np.mod(np.sum(list_v[:,J],axis=0).T,2)!=0).any():
        return False
    A = list_v.T
    index_start_I = 0
    for i in I:
        non_zeros = np.where(A[i]==1)[0]
        if len(non_zeros)==1:
            if non_zeros[0]>index_start_I:
                index_start_I = non_zeros[0]
    index_start_J = 0
    for j in J:
        non_zeros = np.where(A[j]==1)[0]
        if non_zeros[0]>index_start_J:
            index_start_J = non_zeros[0]
    vect = np.zeros(nbr_generators-index_start_I)
    base = 2*np.ones(nbr_generators-index_start_I)
    x = np.ones(nbr_generators)
    for k in range(2**(nbr_generators-index_start_I)):
        x[index_start_I:] = vect
        Ax = np.mod(np.dot(A,x),2)
        if np.count_nonzero(Ax)>G_vector[n-1]:
            give_next_vect(vect,base)
            continue
        R = np.dot(M,Ax)
        if (R==4).any():
            give_next_vect(vect,base)
            continue
        list_results.append(np_K_max_facets[Ax==1].tolist())



# results_path = 'final_results/CSPLS_10_14'
# list_10_14_seeds = [json.loads(facets_bytes) for facets_bytes in read_file(results_path)]
list_results = []
# for K_facets in list_10_14_seeds:
#     # print(K_facets)
#     K = sc.PureSimplicialComplex(K_facets)
#     K_sus = sc.suspension(K)
#     list_IDCM_10_14 = sc.IDCM_Garrison_Scott(K)
#     for IDCM_10_14 in list_IDCM_10_14:
#         IDCM_11_15 = []
#         for x in range(1,16):
#             if x not in IDCM_10_14:
#                 IDCM_11_15.append(x)
#                 break
#         for x in IDCM_10_14:
#             IDCM_11_15.append(x)
#         K_max_facets = sc.find_facets_compatible_with_lambda(IDCM_11_15,m,n)
#         new_f(K_max_facets,K_sus.facets_bin,list_results)
#         print(len(list_results))
nbr = 0

#need to give an order to have the J appear only once
def Try_J(N,K,J,J_max,IDCM_results,starting_index):
    nbr_IDCM = puzalg.puzzle_algo_V2(K,J,True)
    if nbr_IDCM >0:
        if N==0:
            print(J,J_max)
            IDCM_results.append(sc.multiple_wedge(K,J))
        else:
            for k in range(starting_index,len(J)):
                if J[k]+1<=J_max[k]:
                    new_J = J.copy()
                    new_J[k]+=1
                    Try_J(N-1,K,new_J,J_max,IDCM_results,k)



list_IDCM_test = []
for n_seed in range(9,1,-1):
    print(n_seed)
    results_path = 'final_results/CSPLS_%d_%d' % (n_seed,n_seed+p)
    list_m_n_seeds = [json.loads(facets_bytes) for facets_bytes in read_file(results_path)]
    for K_facets in list_m_n_seeds:
        K = sc.PureSimplicialComplex(K_facets)
        J_max = puzalg.is_prediagram_disconnected(K,True)
        if np.sum(J_max)-1 <10:
            continue
        else:
            Try_J(10-n_seed,K,np.zeros(n_seed+p,dtype=int),J_max,list_IDCM_test,0)

list_isom =[list_IDCM_test[0]]
for K in list_IDCM_test[1:]:
    is_isom = False
    for L in list_isom:
        if sc.are_isom(K,L):
            is_isom = True
            break
    if not is_isom:
        list_isom.append(K)
    else:
        del K

data_to_text = []
for K in list_isom:
    print(K.facets_bin)
    data_to_text.append(K.facets_bin)
# text(data_to_text, 'final_results/CnSPLS_10_14')
# N=len(list_isom)
# for k in range(N):
#     K = list_isom[k]
#     list_IDCM_10_14 = sc.IDCM_Garrison_Scott(K)
#     print((k/N)*100,'%',"Nbr of IDCM",len(list_IDCM_10_14))
#     K_sus = sc.suspension(K)
#     if len(list_IDCM_10_14)<100:
#         for IDCM_10_14 in list_IDCM_10_14:
#             IDCM_11_15 = []
#             for x in range(1,16):
#                 if x not in IDCM_10_14:
#                     IDCM_11_15.append(x)
#                     break
#             for x in IDCM_10_14:
#                 IDCM_11_15.append(x)
#             K_max_facets = sc.find_facets_compatible_with_lambda(IDCM_11_15,m,n)
#             new_f(K_max_facets,K_sus.facets_bin,list_results)
#             print(len(list_results))
#     else:
#         print(K.facets_bin)

# text(list_results,raw_results_PATH)
#
# print("Finished",len(list_IDCM_test))
