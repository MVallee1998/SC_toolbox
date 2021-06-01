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
    while index < vect.size - 1 and vect[index] == 0:
        index += 1
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
    ninja = 0
    for IDCM in list_IDCM:
        CM = np.zeros((n, m))
        CM[:, :n] = np.eye(n)
        CM[:, n:m] = IDCM[:n, :]
        list_CM.append(CM.copy())
        print(ninja)
        print(CM)
        ninja+=1
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
                print(i,j,phi)
                vect_of_non_zeros = np.zeros(n)
                vect_of_non_zeros[non_zeros] = 1
                for k in range(m):
                    if (list_CM[i][:, k] == vect_of_non_zeros).all():
                        prediagram_maps[i, j] = phi.copy()
                        prediagram_colors[i, j, k] = 1
                        prediagram_maps[j, i] = phi.copy()
                        prediagram_colors[j, i, k] = 1
    def all_edges_exist(list_generators, list_vertex_of_generator, pos_generator, current_IDCM,
                        current_IDCM_index):
        if pos_generator == len(list_generators):
            return True
        exit_condition = all_edges_exist(list_generators, list_vertex_of_generator, pos_generator + 1,
                                         current_IDCM, current_IDCM_index)
        for index_generator in range(pos_generator, len(list_generators)):
            new_IDCM = np.mod(np.dot(list_generators[index_generator], current_IDCM), 2)
            index_new_IDCM = -1
            for index_CM in range(len(list_CM)):
                if (np.mod(np.dot(list_CM[index_CM], new_IDCM), 2) == 0).all():
                    index_new_IDCM = index_CM
                    break
            if index_new_IDCM < 0:
                print("hello")
                return False
            if prediagram_colors[current_IDCM_index, index_new_IDCM, list_vertex_of_generator[index_generator]] == 0:
                return False
            phi = prediagram_maps[current_IDCM_index,index_new_IDCM]
            indexes_phi = np.flatnonzero(list_CM[current_IDCM_index][:, list_vertex_of_generator[index_generator]])
            generator = np.eye(m)
            for index_pos_phi in indexes_phi:
                generator[index_pos_phi, n:] = phi
            print(current_IDCM_index,index_new_IDCM,list_vertex_of_generator[index_generator])
            print(generator)
            # if (generator!=list_generators[index_generator]).any():
            #     return False
            exit_condition = exit_condition and all_edges_exist(list_generators, list_vertex_of_generator,
                                                                index_generator + 1, new_IDCM, index_new_IDCM)
            if not exit_condition:
                return False
        return exit_condition


    list_conn_comp = []



    for k in range(m):
        list_conn_comp.append([])
        for v in np.unique(prediagram_colors[:, :, k], axis=0):
            list_conn_comp[k].append(np.flatnonzero(v).tolist())
        print(list_conn_comp[k])
    index_wedged_vertices = np.flatnonzero(J_array)
    if index_wedged_vertices.size == 0:
        return list_IDCM_bin
    number_of_realizable_puzzle = 0
    # we choose one IDCM lambda0 for the vertex i_0
    i_0 = index_wedged_vertices[0]
    for i in range(number_IDCM):
        multi_list_of_cases = []
        list_conn_comp = []
        # we then need to enumerate every possible choice of IDCM as a neighbour of lambda0
        for k in index_wedged_vertices:
            conn_comp_k = np.unique(prediagram_colors[:, :, k], axis=0)
            k_conn_comp_of_i = np.flatnonzero(conn_comp_k[conn_comp_k[:, i] == 1])
            list_conn_comp.append(k_conn_comp_of_i)
            cases = np.zeros((k_conn_comp_of_i.size ** (J_array[k]), J_array[k]), dtype=int)
            for l in range(1, cases.shape[0]):
                cases[l] = cases[l - 1].copy()
                give_next_vect(cases[l], k_conn_comp_of_i.size * np.ones(J_array[k]))
            multi_list_of_cases.append(cases)
        number_of_cases = 1
        for v in multi_list_of_cases:
            number_of_cases *= v.shape[0]
        array_of_all_cases = np.zeros((number_of_cases, len(multi_list_of_cases)), dtype=int)
        base_table = [v.shape[0] for v in multi_list_of_cases]
        for k in range(1, number_of_cases):
            array_of_all_cases[k] = array_of_all_cases[k - 1].copy()
            give_next_vect(array_of_all_cases[k], base_table)
        # We have an array with all the indexes of the all the cases
        for case in array_of_all_cases:
            explicit_case = []
            for a in range(case.size):
                explicit_case.append(list_conn_comp[a][multi_list_of_cases[a][case[a]]].tolist())
            # here we list all the generators of the subgroup we chose
            # print(i,explicit_case)
            list_generators = []
            list_vertex_of_generator = []
            for index_list_connected_IDCM in range(len(explicit_case)):
                for index_IDCM in explicit_case[index_list_connected_IDCM]:
                    if i != index_IDCM:
                        phi = prediagram_maps[i, index_IDCM]
                        indexes_phi = np.flatnonzero(list_CM[i][:, index_wedged_vertices[index_list_connected_IDCM]])
                        generator = np.eye(m)
                        for index_pos_phi in indexes_phi:
                            generator[index_pos_phi, n:] = phi
                        already_here = False
                        for other_generator in list_generators:
                            if (generator == other_generator).all():
                                already_here = True
                                break
                        if not already_here:
                            list_generators.append(generator.copy())
                            list_vertex_of_generator.append(index_wedged_vertices[index_list_connected_IDCM])
            # print(list_generators)
            print(i,explicit_case)
            # we then have to try if the subgroup generated by the generators produce only IDCM when acting on lambda0
            number_generators = len(list_generators)
            if all_edges_exist(list_generators, list_vertex_of_generator, 0,
                               list_IDCM[i].copy(), i):
                number_of_realizable_puzzle += 1

            # array_cases_generators = np.zeros((2 ** number_generators, number_generators))
            # for k in range(1, 2 ** number_generators):
            #     array_cases_generators[k] = array_cases_generators[k - 1].copy()
            #     give_next_vect(array_cases_generators[k], 2 * np.ones(number_generators))
            # puzzle_is_realizable = True
            # # print(i,list_generators)
            # for case_generator in array_cases_generators:
            #     if np.count_nonzero(case_generator)>1:
            #         new_IDCM = list_IDCM[i].copy()
            #         for generator_index in np.flatnonzero(case_generator):
            #             new_IDCM = np.mod(np.dot(list_generators[generator_index], new_IDCM), 2)
            #         # print(i,explicit_case,new_IDCM)
            #         IDCM_is_non_singular = False
            #         for index_CM in range(len(list_CM)):
            #             CM = list_CM[index_CM]
            #             if (np.mod(np.dot(CM, new_IDCM), 2) == 0).all():
            #                 IDCM_is_non_singular = True
            #                 break
            #         if not IDCM_is_non_singular:
            #             puzzle_is_realizable = False
            #             break
            # if puzzle_is_realizable:
            #     number_of_realizable_puzzle += 1
            #     if explicit_case[0][0] != explicit_case[1][0] and i!= explicit_case[0][0] and i!= explicit_case[1][0]:
            #         print("hello")
            #         print(i,list_CM[i])
            #         print(explicit_case[0][0],list_CM[explicit_case[0][0]])
            #         print(explicit_case[1][0],list_CM[explicit_case[1][0]])
    print(number_of_realizable_puzzle)


K = sc.PureSimplicialComplex([[1, 2], [1, 6], [2, 3], [3, 4], [4, 5], [5, 6]])
print("nbr1")
puzzle_algo(K, [1, 0, 0, 0, 0, 1])
# print("nbr2")
# puzzle_algo(K, [1, 1, 0, 0, 0, 0])
wedged_K1 = sc.multiple_wedge(K, [1, 1 ,0, 0, 0, 0])
print(len(sc.Garrison_Scott(wedged_K1)))

# for n in range(2, 3):
#     m = n + 4
#     results = read_file('final_results_BAK/PLS_%d_%d' % (m, n))
#     start = timeit.default_timer()
#     for K_byte in results:
#         K = sc.PureSimplicialComplex(json.loads(K_byte))
#         puzzle_algo(K, [1,0,0,1,0,0])
#         print(len(sc.Garrison_Scott(sc.multiple_wedge(K,[1,0,0,1,0,0]))))
