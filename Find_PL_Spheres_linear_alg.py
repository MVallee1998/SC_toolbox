import linear_alg_method as lam
import timeit
import numpy as np
import cupy as cp
import SimplicialComplex as sc
import numba as nb
from itertools import combinations
import scipy.sparse
import sys
np.set_printoptions(threshold=sys.maxsize)

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]

np_arrange = np.arange(0,256)
np_arrange_odd = 2*np.arange(0,127) + 1
m = 11
n = 7
number_steps = 6
raw_results_PATH = 'raw_results/all_PLS_%d_%d' % (m, n)

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

def Gauss(M):
    N = M.copy()
    a,b = M.shape
    i = 0
    j = 0
    while i<a and j<b:
        if N[i,j]==0:
            i_0 = -1
            for i_iter in range(i+1,a):
                if N[i_iter,j] == 1:
                    i_0 = i_iter
                    break
            if i_0 == -1:
                j+=1
                continue
            N[[i_0, i]] = N[[i, i_0]]
        for i_iter in range(0,i):
            N[i_iter] = (N[i_iter] + N[i]*N[i_iter,j]) % 2
        for i_iter in range(i + 1, a):
            if N[i_iter, j] == 1:
                N[i_iter] = (N[i_iter] + N[i]) % 2
        j += 1
        i += 1
    return(N)

def reduce_wrt_columns(M,array_columns,starting_row):
    N = M.copy()
    a,c = M.shape
    b = array_columns.size
    i = starting_row
    j = 0
    while i<a and j<b:
        if N[i,array_columns[j]]==0:
            i_0 = -1
            for i_iter in range(i+1,a):
                if N[i_iter,array_columns[j]] == 1:
                    i_0 = i_iter
                    break
            if i_0 == -1:
                j+=1
                continue
            N[[i_0, i]] = N[[i, i_0]]
        for i_iter in range(i):
            N[i_iter] = (N[i_iter] + N[i]*N[i_iter,array_columns[j]]) % 2
        for i_iter in range(i + 1, a):
            N[i_iter] = (N[i_iter] + N[i]*N[i_iter, array_columns[j]]) % 2
        j += 1
        i += 1
    return(N)


def rank(M):
    return(np.count_nonzero(np.sum(Gauss(M),axis=1)))

@nb.njit
def give_next_vect(vect, base):
    index = 0
    vect[index] = (vect[index] + 1) % base[index]
    while index < vect.size - 1 and vect[index] == 0:
        index += 1
        vect[index] = (vect[index] + 1) % base[index]
    if index == vect.size -1:
        return False
    return True


def f(char_funct):
    start = timeit.default_timer()
    facets = sc.find_facets_compatible_with_lambda(char_funct,m,n)
    M = lam.construct_matrix(facets)
    M_cp = cp.asarray(M)
    list_v = lam.find_kernel(M)
    nbr_results = list_v.shape[0]
    print(nbr_results)
    print(M.shape)
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
        both_condition = cp.logical_and(verifying_G_theorem, having_every_closed_ridges)
        if cp.sum(both_condition):
            good_conditions = cp.flatnonzero(both_condition)
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


def new_f(facets):
    start = timeit.default_timer()
    M = lam.construct_matrix(facets)
    M_cp = cp.asarray(M)
    M_sparse = scipy.sparse.csr_matrix(M)
    list_v = lam.find_kernel(M)
    list_v_new = reduce_wrt_columns(list_v,np.array([0]),0)
    nbr_results = list_v.shape[0]
    #The idea is to reorganize the generators so some subset of them cannot be added together
    list_not_together = M[np.sum(M,axis = 1)==5]
    sum_of_not_together = np.zeros(M.shape[1]) # this array represents which MF have been used already
    sum_of_not_together+=list_v_new[0,:]
    list_distinct_not_together = []
    starting_row = 1
    list_gener_not_together = []
    find_new_one= True
    while find_new_one:
        find_new_one = False
        for not_together in list_not_together:
            index_not_together = np.flatnonzero(not_together)
            if np.sum(np.multiply(not_together,sum_of_not_together))==0 and rank(list_v_new[starting_row:,index_not_together]) == index_not_together.size -1:
                sum_of_not_together+=not_together
                list_v_new = reduce_wrt_columns(list_v_new,index_not_together,starting_row)
                list_gener_not_together.append(list(range(starting_row,starting_row+index_not_together.size-1)))
                starting_row+=index_not_together.size-1
                list_distinct_not_together.append(not_together.copy())
                find_new_one = True
                break
    list_not_together = M[np.sum(M,axis = 1)==4]
    find_new_one= True
    while find_new_one:
        find_new_one = False
        for not_together in list_not_together:
            index_not_together = np.flatnonzero(not_together)
            if np.sum(np.multiply(not_together,sum_of_not_together))==0 and rank(list_v_new[starting_row:,index_not_together]) == index_not_together.size -1:
                sum_of_not_together+=not_together
                list_v_new = reduce_wrt_columns(list_v_new,index_not_together,starting_row)
                list_gener_not_together.append(list(range(starting_row,starting_row+index_not_together.size-1)))
                starting_row+=index_not_together.size-1
                list_distinct_not_together.append(not_together.copy())
                find_new_one = True
                break
    for k in range(starting_row,nbr_results):
        list_gener_not_together.append([k])
    #here we will create the lists where we store how to build the linear sums
    list_to_pick_lin_comb = []
    array_number_lines = np.zeros(len(list_gener_not_together),dtype=np.uint64)
    number_cases = 1
    for index in range(len(list_gener_not_together)):
        not_together = list_gener_not_together[index]
        if len(not_together)==1:
            nbr_lines = 2
        elif len(not_together)==3:
            nbr_lines = 7
        else:
            nbr_lines = 11
        array_number_lines[index] = nbr_lines
        list_to_pick_lin_comb.append(np.zeros((nbr_lines,nbr_results)))
        current_line=1
        for k in range(1,3):
            for iter_combi in combinations(not_together,k):
                list_to_pick_lin_comb[-1][current_line,list(iter_combi)] = 1
                current_line+=1
        number_cases*=nbr_lines
    base_vect_to_mult_array = np.zeros((np.prod(array_number_lines[:number_steps]),nbr_results))
    base_vect_to_mult_array[:,0] = 1
    print(nbr_results,np.format_float_scientific(np.prod(array_number_lines)))
    vect = np.zeros(number_steps,dtype=int)
    for k in range(1,np.prod(array_number_lines[:number_steps])):
        give_next_vect(vect,array_number_lines[:number_steps])
        for l in range(number_steps):
            base_vect_to_mult_array[k] += list_to_pick_lin_comb[l][vect[l],:]
    np_facets = cp.array(facets)
    A = cp.asarray(np.transpose(list_v_new))
    results = []
    vect = np.zeros(len(array_number_lines)-number_steps,dtype=int)
    keep_going = True

    #this is the main loop
    while keep_going:
        vect_to_mult_array = base_vect_to_mult_array.copy()
        for l in range(number_steps,number_steps+vect.size):
            vect_to_mult_array += list_to_pick_lin_comb[l][vect[l-number_steps]]
        candidate_array, prod = get_product(M_cp,A, cp.asarray(vect_to_mult_array.T))
        verifying_G_theorem = cp.sum(candidate_array, axis=0) <= G_vector[n - 1]
        having_every_closed_ridges = cp.logical_not((prod >= 4).any(axis=0))
        both_condition = cp.logical_and(verifying_G_theorem, having_every_closed_ridges)
        good_conditions = cp.flatnonzero(both_condition)
        good_candidates = candidate_array[:,good_conditions].T
        for good_candidate in good_candidates:
            good_candidate_facets = np_facets[good_candidate==1]
            good_candidate_facets_list = good_candidate_facets.tolist()
            results.append(good_candidate_facets_list)
            # K = sc.PureSimplicialComplex(good_candidate_facets_list)
            # text(good_candidate_facets_list,raw_results_PATH)
        keep_going = give_next_vect(vect,array_number_lines[number_steps:])
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    # print("number of results", len(results))
    return results


# if __name__ == '__main__':
#     list_char_funct = sc.enumerate_char_funct_orbits(n, m)
#     with Pool(processes=6) as pool:
#         big_result = pool.imap(f, list_char_funct)
#         for results in big_result:
#             text(results,raw_results_PATH)

list_char_funct = sc.enumerate_char_funct_orbits(n, m)
for char_funct in list_char_funct[:1]:
    facets = sc.find_facets_compatible_with_lambda(char_funct,m,n)
    results = new_f(facets)


# for n in range(2,8):
#     m=n+4
#     list_char_funct = sc.enumerate_char_funct_orbits(n, m)
#     global_start = timeit.default_timer()
#     if n>4:
#         number_steps = 4
#     for char_funct in list_char_funct[:1]:
#         facets = sc.find_facets_compatible_with_lambda(char_funct,m,n)
#         results = new_f(facets)
#         # text(results,raw_results_PATH)
#     global_end = timeit.default_timer()
#     print((n,m), global_end - global_start)


# for n in range(2,5):
#     m=n+4
#     MFset = []
#     print(n,m)
#     for MF in combinations(range(1, m + 1), n):
#         MFset.append(sc.face_to_binary(MF,m))
#     results = new_f(MFset)
# global_end = timeit.default_timer()
# text(results,raw_results_PATH)
# print("Total time spent", global_end - global_start)