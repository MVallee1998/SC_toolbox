import numpy as np
import timeit
import SimplicialComplex as sc
import json

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




def puzzle_algo_V2(K,J):
    J_array = np.array(J)
    p = K.Pic
    n = K.n
    m = K.m
    list_IDCM_bin = sc.IDCM_Garrison_Scott(K)
    # print("number of IDCM: ", len(list_IDCM_bin))
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
                vect_of_non_zeros = np.zeros(n)
                vect_of_non_zeros[non_zeros] = 1
                for k in range(m):
                    if (list_CM[i][:, k] == vect_of_non_zeros).all():
                        prediagram_maps[i, j] = phi.copy()
                        prediagram_colors[i, j, k] = 1
                        prediagram_maps[j, i] = phi.copy()
                        prediagram_colors[j, i, k] = 1
    nbr_vertices = np.prod(J_array+1)
    array_vertices = np.zeros((nbr_vertices,m))
    for k in range(1,nbr_vertices):
        array_vertices[k] = array_vertices[k-1].copy()
        give_next_vect(array_vertices[k],J_array+1)
    array_depth = (np.count_nonzero(array_vertices,axis = 1))
    list_indexes_depth = []
    for depth in range(np.count_nonzero(J_array)+1):
        list_indexes_depth.append(np.flatnonzero(array_depth==depth).tolist())
    def construct_puzzle(depth, position, current_puzzle,list_puzzles):
        if depth > np.count_nonzero(J_array):
            list_puzzles.append(current_puzzle)
        elif position == len(list_indexes_depth[depth]):
            construct_puzzle(depth+1, 0, current_puzzle,list_puzzles)
        else:
            current_vertex = array_vertices[list_indexes_depth[depth][position]]
            list_neighbours = []
            # I find all the vertices adjacent to the current vertex whi has already their images fixed
            for lower_depth in range(depth+1):
                for index_other_vertex in list_indexes_depth[lower_depth]:
                    if lower_depth<depth or index_other_vertex < list_indexes_depth[depth][position]:
                        other_vertex = array_vertices[index_other_vertex]
                        position_of_edge = np.flatnonzero(current_vertex - other_vertex)
                        if position_of_edge.size == 1:
                            list_neighbours.append((index_other_vertex,position_of_edge[0]))
            marker_of_possible_IDCM = np.ones(number_IDCM)
            for index_data_neighbour in range(len(list_neighbours)):
                data_neighbour = list_neighbours[index_data_neighbour]
                index_neighbour , p = data_neighbour
                neighbour_IDCM = current_puzzle[index_neighbour]
                for data_neighbour2 in list_neighbours[index_data_neighbour+1:]:
                    index_neighbour2, q = data_neighbour2
                    neighbour_IDCM2 = current_puzzle[index_neighbour2]
                    if p!=q and (list_CM[neighbour_IDCM][:,p] == list_CM[neighbour_IDCM][:,q]).all() and  (list_CM[neighbour_IDCM2][:,p] == list_CM[neighbour_IDCM2][:,q]).all() and  (list_CM[neighbour_IDCM][:,p] == list_CM[neighbour_IDCM2][:,q]).all():
                        opposite_corner_vertex = array_vertices[index_neighbour].copy()
                        opposite_corner_vertex[q] = array_vertices[index_neighbour2][q]
                        index_opposite_corner = -1
                        for k in range(nbr_vertices):
                            if (opposite_corner_vertex==array_vertices[k]).all():
                                index_opposite_corner = k
                                break
                        opposite_corner_IDCM = current_puzzle[index_opposite_corner]
                        mandatory_IDCM = np.mod(list_IDCM[opposite_corner_IDCM]+list_IDCM[neighbour_IDCM]+ list_IDCM[neighbour_IDCM2],2)
                        index_mandatory_IDCM = -1
                        for index_IDCM in range(len(list_IDCM)):
                            if (list_IDCM[index_IDCM]==mandatory_IDCM).all():
                                index_mandatory_IDCM = index_IDCM
                                break
                        marker_mandatory_IDCM = np.zeros(number_IDCM)
                        if index_mandatory_IDCM>=0:
                            marker_mandatory_IDCM[index_mandatory_IDCM] = 1
                        marker_of_possible_IDCM = np.multiply(marker_of_possible_IDCM,marker_mandatory_IDCM)
                        if np.count_nonzero(marker_of_possible_IDCM) == 0:
                            break
                marker_of_possible_IDCM = np.multiply(marker_of_possible_IDCM, prediagram_colors[:,neighbour_IDCM,p])
                if np.count_nonzero(marker_of_possible_IDCM) == 0:
                    break
            indexes_possible_IDCM = np.flatnonzero(marker_of_possible_IDCM)
            if indexes_possible_IDCM.size != 0:
                for index_IDCM in indexes_possible_IDCM:
                    new_puzzle = current_puzzle.copy()
                    new_puzzle[list_indexes_depth[depth][position]] = index_IDCM
                    construct_puzzle(depth,position+1,new_puzzle,list_puzzles)
    original_puzzle = -1*np.ones(nbr_vertices,dtype = int)
    final_list_puzzles = []
    construct_puzzle(0,0,original_puzzle,final_list_puzzles)
    if len(final_list_puzzles)!= len(list_IDCM):
        K.compute_MNF_set()
        K.MNF_bin_to_MNF()
        print("hello",K.MNF_set,len(final_list_puzzles), len(list_IDCM))
    # return len(final_list_puzzles)




# K = sc.PureSimplicialComplex([[1, 2], [1, 6], [2, 3], [3, 4], [4, 5],[5,6]])
# wedged_K1 = sc.multiple_wedge(K, [6,1 ,0, 0, 0,0])
# print(wedged_K1.facets)
# start = timeit.default_timer()
# print(len(sc.Garrison_Scott(wedged_K1)))
# stop = timeit.default_timer()
# print("Time spent for GS: ", stop - start)
# start = timeit.default_timer()
# puzzle_algo_V2(K, [6, 1, 0, 0, 0,0])
# stop = timeit.default_timer()
# print("Time spent for puzzle: ", stop - start)


for n in range(2, 9):
    m = n + 4
    results = read_file('final_results/PLS_%d_%d' % (m, n))
    start = timeit.default_timer()
    for i in range(len(results)):
        K_byte = results[i]
        K = sc.PureSimplicialComplex(json.loads(K_byte))
        # K.compute_MNF_set()
        # K.MNF_bin_to_MNF()
        J = [0]*K.m
        J[0]=1
        puzzle_algo_V2(K,J)
        # print(i/len(results)*100,"%")
    stop = timeit.default_timer()
    print("(",n,",",m,")","Puzzle mean",(stop-start)/len(results))
    # start = timeit.default_timer()
    # for i in range(min(len(results),100)):
    #     K_byte = results[i]
    #     K = sc.PureSimplicialComplex(json.loads(K_byte))
    #     J = [1]*K.m
    #     print(i/min(len(results),100)*100,"%")
    #     sc.Garrison_Scott(sc.multiple_wedge(K,J))
    # stop = timeit.default_timer()
    # print("(",n,",",m,")","GS mean",(stop-start)/min(len(results),100))



