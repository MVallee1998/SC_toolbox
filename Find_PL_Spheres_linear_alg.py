import linear_alg_method as lam
import timeit
import numpy as np
import numpy as cp
import SimplicialComplex as sc
import numba as nb
from itertools import combinations
import sys

np.set_printoptions(threshold=sys.maxsize)

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]

np_arrange = np.arange(0, 256)
np_arrange_odd = 2 * np.arange(0, 127) + 1
m = 8
n = 4
number_steps = 1

raw_results_PATH = 'raw_results/PLS_%d_%d' % (m, n)


def text(results, path):
    t = open(path, mode='a', encoding='utf-8')
    for K in results:
        t.write(str(K) + '\n')
    t.close()


@nb.njit
def get_product(M, A, vect_to_mult_array):
    candidate_array = A.dot(vect_to_mult_array) % 2
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
    while vect[index] == 0 and index < vect.size-1:
        index += 1
        vect[index] = (vect[index] + 1) % base[index]
    if index == vect.size-1:
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


def new_f(facets):
    start = timeit.default_timer()
    M = lam.construct_matrix(facets)
    M_cp = cp.asarray(M)
    list_v = lam.find_kernel(M)
    reduce_wrt_columns(list_v, np.array([0]), 0)
    nbr_results = list_v.shape[0]

    # The idea is to reorganize the generators so some subset of them cannot be added together

    sum_of_not_together = np.zeros(M.shape[1])  # this array represents which MF have been used already
    sum_of_not_together += list_v[0, :]
    list_distinct_not_together = []
    starting_row = 1
    list_gener_not_together = []
    list_not_together = M[np.sum(M, axis=1) == 6]
    starting_row = partition_generators(list_v, starting_row, list_not_together, list_distinct_not_together,
                                        sum_of_not_together,
                                        list_gener_not_together)
    list_not_together = M[np.sum(M, axis=1) == 5]
    starting_row = partition_generators(list_v, starting_row, list_not_together, list_distinct_not_together,
                                        sum_of_not_together,
                                        list_gener_not_together)
    list_not_together = M[np.sum(M, axis=1) == 4]
    starting_row = partition_generators(list_v, starting_row, list_not_together, list_distinct_not_together,
                                        sum_of_not_together,
                                        list_gener_not_together)

    for k in range(starting_row, nbr_results):
        list_gener_not_together.append([k])
    # here we will create the lists where we store how to build the linear sums
    list_to_pick_lin_comb = []
    array_number_lines = np.zeros(len(list_gener_not_together), dtype=np.uint64)
    number_cases = 1
    for index in range(len(list_gener_not_together)):
        not_together = list_gener_not_together[index]
        if len(not_together) == 1:
            nbr_lines = 2
        elif len(not_together) == 3:
            nbr_lines = 7
        elif len(not_together) == 4:
            nbr_lines = 11
        else:
            nbr_lines = 16
        array_number_lines[index] = nbr_lines
        list_to_pick_lin_comb.append(np.zeros((nbr_lines, nbr_results)))
        current_line = 1
        for k in range(1, 3):
            for iter_combi in combinations(not_together, k):
                list_to_pick_lin_comb[-1][current_line, list(iter_combi)] = 1
                current_line += 1
        number_cases *= nbr_lines
    base_vect_to_mult_array = np.zeros((np.prod(array_number_lines[:number_steps]), nbr_results))
    base_vect_to_mult_array[:, 0] = 1
    print(nbr_results, array_number_lines, np.format_float_scientific(np.prod(array_number_lines)))
    vect = np.zeros(number_steps, dtype=int)
    for k in range(1, np.prod(array_number_lines[:number_steps])):
        give_next_vect(vect, array_number_lines[:number_steps])
        for l in range(number_steps):
            base_vect_to_mult_array[k] += list_to_pick_lin_comb[l][vect[l], :]
    np_facets = cp.array(facets)
    A = cp.asarray(np.transpose(list_v))

    results = []
    vect = np.zeros(len(array_number_lines) - number_steps, dtype=int)
    keep_going = True
    # this is the main loop
    while keep_going:
        vect_to_mult_array = base_vect_to_mult_array.copy()
        for l in range(number_steps, number_steps + vect.size):
            vect_to_mult_array += list_to_pick_lin_comb[l][vect[l - number_steps]]
        candidate_array, prod = get_product(M_cp, A, cp.asarray(vect_to_mult_array.T))
        # verifying_G_theorem = cp.sum(candidate_array, axis=0) >-10
        verifying_G_theorem = cp.sum(candidate_array, axis=0) <= G_vector[n - 1]
        having_every_closed_ridges = cp.logical_not((prod >= 4).any(axis=0))
        both_condition = cp.logical_and(verifying_G_theorem, having_every_closed_ridges)
        good_conditions = cp.flatnonzero(both_condition)
        good_candidates = candidate_array[:, good_conditions].T
        for good_candidate in good_candidates:
            good_candidate_facets = np_facets[good_candidate == 1]
            good_candidate_facets_list = good_candidate_facets.tolist()
            results.append(good_candidate_facets_list)
            # K = sc.PureSimplicialComplex(good_candidate_facets_list)
            # text(good_candidate_facets_list,raw_results_PATH)
        keep_going = give_next_vect(vect, array_number_lines[number_steps:])
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    return results


# if __name__ == '__main__':
#     list_char_funct = sc.enumerate_char_funct_orbits(n, m)
#     with Pool(processes=6) as pool:
#         big_result = pool.imap(f, list_char_funct)
#         for results in big_result:
#             text(results,raw_results_PATH)

# list_char_funct = sc.enumerate_char_funct_orbits(n, m)
#
# for char_funct in list_char_funct:
#     facets = sc.find_facets_compatible_with_lambda(char_funct, m, n)
#     results = new_f(facets)
#     # text(results, raw_results_PATH)

for n in range(2,8):
    m=n+4
    list_char_funct = sc.enumerate_char_funct_orbits(n, m)
    global_start = timeit.default_timer()
    if n>4:
        number_steps = 3
    for char_funct in list_char_funct:
        facets = sc.find_facets_compatible_with_lambda(char_funct,m,n)
        results = new_f(facets)
        raw_results_PATH = 'raw_results/PLS_%d_%d' % (m, n)
        text(results,raw_results_PATH)
    global_end = timeit.default_timer()
    print((n,m), global_end - global_start)


# for n in range(3,4):
#     m=n+3
#     MFset = []
#     print(n,m)
#     for MF in combinations(range(1, m + 1), n):
#         MFset.append(sc.face_to_binary(MF,m))
#     results = new_f(MFset)
#     raw_results_PATH = 'raw_results/all_PLS_%d_%d' % (m, n)
#     text(results, raw_results_PATH)

# print("Total time spent", global_end - global_start)

print("Finished")