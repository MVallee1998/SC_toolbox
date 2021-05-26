import json
import timeit

import numpy as np

import SimplicialComplex as sc

list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1] * 2)


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

def give_next_vect(vect, base):
    index = 0
    vect[index] = (vect[index] + 1) % base[index]
    while index<vect.size-1 and vect[index]==0:
        index+=1
        vect[index] = (vect[index] + 1) % base[index]


def puzzle_algo(K, J):
    J_array = np.array(J)
    p = K.Pic
    n = K.n
    m = K.m
    list_IDCM_bin = sc.IDCM_Garrison_Scott(K)
    list_IDCM = []
    for IDCM_bin in list_IDCM_bin:
        IDCM = np.zeros((m, p))
        for i in range(m):
            for j in range(p):
                if list_2_pow[j] | IDCM_bin[i] == IDCM_bin[i]:
                    IDCM[i, j] = 1
        list_IDCM.append(IDCM.copy())
    list_CM = []
    for IDCM in list_IDCM:
        CM = np.zeros((n, m))
        CM[:, :n] = np.eye(n)
        CM[:, n:m] = IDCM[:n, :]
        list_CM.append(CM.copy())
    number_IDCM = len(list_IDCM)
    # constructing the prediagram
    prediagram_colors = np.zeros((number_IDCM, number_IDCM, m))
    prediagram_maps = np.zeros((number_IDCM, number_IDCM, p))
    for k in range(number_IDCM):
        prediagram_colors[k, k] = np.ones(m)
    for i in range(number_IDCM):
        for j in range(i + 1, number_IDCM):
            sum_IDCM = np.mod(list_IDCM[i] + list_IDCM[j], 2)
            non_zeros = np.flatnonzero((sum_IDCM == 1).any(axis=1))
            unique_maps = np.unique(sum_IDCM[non_zeros], axis=0)
            if (unique_maps.shape[0]) == 1:
                phi = unique_maps[0]
                vect_of_non_zeros = np.zeros(n)
                vect_of_non_zeros[non_zeros] = 1
                for k in range(m):
                    if (list_CM[i][:, k] == vect_of_non_zeros).all():
                        if k < n:
                            prediagram_maps[i, j] = phi.copy()
                            prediagram_colors[i, j, k] = 1
                            prediagram_maps[j, i] = phi.copy()
                            prediagram_colors[j, i, k] = 1
                        else:
                            if phi[k - n] == 0:
                                prediagram_maps[i, j] = phi.copy()
                                prediagram_colors[i, j, k] = 1
                                prediagram_maps[j, i] = phi.copy()
                                prediagram_colors[j, i, k] = 1
    list_conn_comp = []
    for k in range(p):
        list_conn_comp.append([])
        for v in np.unique(prediagram_colors[:, :, k], axis=0):
            list_conn_comp[k].append(np.flatnonzero(v).tolist())
    index_wedged_vertices = np.flatnonzero(J_array)
    if index_wedged_vertices.size == 0:
        return list_IDCM_bin
    #we choose one IDCM lambda0 for the vertex i_0
    for i in range(number_IDCM):
        multi_list_of_cases = []
        list_conn_comp = []
        #we then need to enumerate every possible choice of IDCM as a neighbour of lambda0
        for k in index_wedged_vertices:
            conn_comp_k = np.unique(prediagram_colors[:,:,k],axis=0)
            k_conn_comp_of_i = np.flatnonzero(conn_comp_k[conn_comp_k[:,i]==1])
            list_conn_comp.append(k_conn_comp_of_i)
            cases = np.zeros((k_conn_comp_of_i.size**(J_array[k]),J_array[k]),dtype=int)
            for l in range(1,cases.shape[0]):
                cases[l] = cases[l-1].copy()
                give_next_vect(cases[l],k_conn_comp_of_i.size*np.ones(J_array[k]))
            multi_list_of_cases.append(cases)
        number_of_cases = 1
        for v in multi_list_of_cases:
            number_of_cases*=v.shape[0]
        print(number_of_cases)
        array_of_all_cases = np.zeros((number_of_cases,len(multi_list_of_cases)),dtype= int)
        base_table = [v.shape[0] for v in multi_list_of_cases]
        for k in range(1,number_of_cases):
            array_of_all_cases[k] = array_of_all_cases[k-1].copy()
            give_next_vect(array_of_all_cases[k],base_table)
        #We have an array with all the indexes of the all the cases
        for case in array_of_all_cases:
            explicit_case = []
            for a in range(case.size):
                explicit_case.append(list_conn_comp[a][multi_list_of_cases[a][case[a]]].tolist())
            list_generators = []
            for index_list_connected_IDCM in range(len(explicit_case)):
                list_generators.append([])
                for index_IDCM in explicit_case[index_list_connected_IDCM]:
                    if i!=index_IDCM:
                        phi = prediagram_maps[i,index_IDCM]
                        indexes_phi = np.flatnonzero(list_CM[i][:,index_wedged_vertices[index_list_connected_IDCM]])

                        generator = np.eye(m)
                        for index_pos_phi in indexes_phi:
                            generator[index_pos_phi,n:] = phi
                        if indexes_phi.size>1:
                            print( list_IDCM[i])
                            print(list_IDCM[index_IDCM])






K = sc.PureSimplicialComplex([[1,2],[1,6],[2,3],[3,4],[4,5],[5,6]])
puzzle_algo(K, [2,0,0,0,2,0])
print(len(sc.Garrison_Scott(sc.multiple_wedge(K,[2,0,0,0,0,0]))))

# for n in range(2, 3):
#     m = n + 4
#     results = read_file('final_results_BAK/PLS_%d_%d' % (m, n))
#     start = timeit.default_timer()
#     for K_byte in results:
#         K = sc.PureSimplicialComplex(json.loads(K_byte))
#         puzzle_algo(K, [1,0,0,1,0,0])
#         print(len(sc.Garrison_Scott(sc.multiple_wedge(K,[1,0,0,1,0,0]))))