import linear_alg_method as lam
import timeit
import numpy as np
import numpy as cp
import SimplicialComplex as sc
import numba as nb
from itertools import combinations
import sys


pow_2 = np.ones(64,dtype=np.uint64)
for k in range(62,-1,-1):
    pow_2[k] = pow_2[k+1]*2
# cfoo(pow_2)
np.set_printoptions(threshold=sys.maxsize)

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]

np_arrange = np.arange(0, 256)
np_arrange_odd = 2 * np.arange(0, 127) + 1
m = 9
n = 5
p=m-n
number_steps = 1

raw_results_PATH = 'raw_results/CSPLS_%d_%d' % (n,m)


def text(results, path):
    t = open(path, mode='a', encoding='utf-8')
    for K in results:
        t.write(str(K) + '\n')
    t.close()


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


def new_f(facets,index_data):
    M = lam.construct_matrix(facets)
    nbr_ridges, nbr_facets = M.shape
    print(M.shape)
    list_v = lam.find_kernel(M.copy())
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
    list_groups = np.ones(len(list_gener_not_together)+1,dtype=np.uint16)
    for index in range(len(list_gener_not_together)):
        not_together = list_gener_not_together[index]
        list_groups[index+1] = len(not_together)
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
    print("dimension of the kernel of A:",nbr_results)
    print("number of vectors to compute: ", np.count_nonzero(array_number_lines==11), " ", np.count_nonzero(array_number_lines==2))
    vect = np.zeros(number_steps, dtype=int)
    return 0
    for k in range(1, np.prod(array_number_lines[:number_steps])):
        give_next_vect(vect, array_number_lines[:number_steps])
        for l in range(number_steps):
            base_vect_to_mult_array[k] += list_to_pick_lin_comb[l][vect[l], :]
    np_facets = cp.array(facets)
    A = cp.asarray(np.transpose(list_v),dtype=np.uint)
    A_to_C = np.dot(A,pow_2[:nbr_results])
    t = open('CUDA_DATA/Data_'+str(n)+'_'+str(m)+'_'+str(index_data)+'.cpp', mode='a', encoding='utf-8')
    t.write('#define N '+str(n))
    t.write('\n')
    t.write('#define NBR_FACETS '+str(nbr_facets))
    t.write('\n')
    t.write('#define NBR_GENERATORS '+str(nbr_results))
    t.write('\n')
    t.write('#define NBR_X0 '+str(np.prod(array_number_lines[2:8])))
    t.write('\n')
    t.write('#define NBR_X1 '+str(np.prod(array_number_lines[8:])))
    t.write('\n')
    # t.write('#define NBR_RIDGES '+str(((nbr_ridges//196)+1)*196))
    # t.write('\n')
    t.write('#define NBR_GROUPS '+ str(len(list_groups)))
    t.write('\n')
    t.write("#define MAX_NBR_FACETS "+str(G_vector[n-1]))
    t.write('\n')
    t.write("int sizeVectX0 = 8;")
    t.write('\n')
    t.write("int sizeVectX1 = "+str(max(0,len(list_groups)-9))+';')
    t.write('\n')
    t.write('unsigned long A[NBR_FACETS]={')
    for k in range(nbr_facets):
        if k<nbr_facets-1:
            t.write(str(A_to_C[k])+',')
        else:
            t.write(str(A_to_C[k])+'};')
    M_to_C = np.ones((n,nbr_facets),np.uint16)
    t.write('\n')
    for i in range(nbr_facets):
        spot = 0
        for j in range(nbr_ridges):
            if M[j,i]==1:
                M_to_C[spot,i] = j
                spot+=1
    t.write('unsigned int M[N][NBR_FACETS]={')
    for i in range(n):
        t.write('{')
        for j in range(nbr_facets):
            if j<nbr_facets:
                t.write(str(M_to_C[i,j])+',')
            else:
                t.write(str(M_to_C[i,j]))
        t.write('},')
    t.write('};')
    t.write('\n')
    t.write('unsigned int F[NBR_FACETS]={')
    for F in np_facets:
        t.write(str(F)+',')
    t.write('};')
    t.write('\n')
    t.write("unsigned int list_groups[NBR_GROUPS] = {")
    for i in range(len(list_groups)):
        if i<len(list_groups-1):
            t.write(str(list_groups[i])+',')
        else:
            t.write(str(list_groups[i]))
    t.write('};')
    t.write('\n')
    t.close()
    print(list_groups,np.sum(list_groups))


MFset = []
for MF in combinations(range(1, m + 1), n):
    MFset.append(sc.face_to_binary(MF, m))
new_f(MFset,0)
print("Finished")
