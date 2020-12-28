from itertools import combinations
import Z2_linear_algebra
import numpy as np
import SimplicialComplex as sc


def enumerate_facets_and_ridges(char_function, n, m):
    Pic = m - n
    cofacets = []
    for cofacet_iter in combinations(range(1, m + 1), Pic):
        sub_array = []
        for index in cofacet_iter:
            sub_array.append(char_function[index - 1])
        if Z2_linear_algebra.Z2Array(Pic, sub_array.copy()).is_invertible():
            cofacets.append(list(cofacet_iter))
    candidate_facets = []
    list_2_pow = [1]
    for k in range(1, m + 1):
        list_2_pow.append(list_2_pow[-1] * 2)
    for cofacet in cofacets:
        candidate_facets.append((list_2_pow[m] - 1) ^ sc.face_to_binary(cofacet, m))
    candidate_facets.sort()
    ridges = []
    for facet in candidate_facets:
        for element in list_2_pow:
            if element | facet == facet:
                ridge = element ^ facet
                if not ridge in ridges:
                    ridges.append(ridge)
    ridges.sort()
    print("Facets and ridges enumerated")
    return candidate_facets, ridges


def construct_graph(char_funct, n, m):
    full_simplex = sc.list_2_pow[m] - 1
    facets, ridges = enumerate_facets_and_ridges(char_funct, n, m)
    G = [(k, []) for k in range(len(facets))]
    for i in range(len(G)):
        for element1 in sc.list_2_pow:
            if element1 | facets[i] == facets[i]:
                ridge = element1 ^ facets[i]
                position_ridge = sc.dichotomie(ridges, ridge)
                G[i][1].append((position_ridge, []))
                other_vertices = full_simplex ^ facets[i]
                for element2 in sc.list_2_pow:
                    if element2 | other_vertices == other_vertices:
                        position_of_connected_facet = sc.dichotomie(facets, element2 | ridge)
                        if position_of_connected_facet >= 0:
                            G[i][1][-1][1].append(position_of_connected_facet)

        # for j in range(len(facets)):
        #     ridge_index = dichotomie(ridges, facets[i] & facets[j])
        #     if ridge_index >= 0:
        #         G[i][2].append((j,ridge_index))
        #         if ridge_index not in G[i][1]:
        #             G[i][1].append(ridge_index)
    print("Graph constructed")
    return facets, ridges, G


def enumerate_cases(list_of_max, k, current_list, results):
    if k == len(list_of_max):
        results.append(current_list)
    else:
        if list_of_max[k] == -1:
            enumerate_cases(list_of_max, k + 1, current_list.copy(), results)
        else:
            enumerate_cases(list_of_max, k + 1, current_list.copy(), results)
            if current_list[k] < list_of_max[k] - 1:
                new_list = current_list.copy()
                new_list[k] += 1
                enumerate_cases(list_of_max, k, new_list, results)


def find_layer(G, facet_layer, queue, current_layer):
    if -1 in facet_layer:
        new_queue = []
        for facet_index in queue:
            for data_ridges in G[facet_index][1]:
                for neighbor_facet_index in data_ridges[1]:
                    if facet_layer[neighbor_facet_index] == -1:
                        facet_layer[neighbor_facet_index] = current_layer
                        new_queue.append(neighbor_facet_index)
        find_layer(G, facet_layer, new_queue, current_layer + 1)


def graph_method2(facets, ridges, G, max_counter, starting_point=None):
    def graph_method_rec(info_ridges, info_facets, facet_layer, queued_facets, counter, nbr_facets, nbr_ridges,
                         result_K, max_counter):
        if counter <= 6:
            if max_counter == 0 and queued_facets == []:
                if 1 not in info_ridges:
                    K = []
                    for k in range(len(info_facets)):
                        if info_facets[k] == 1:
                            K.append(facets[k])
                    K.sort()
                    # print(K,counter)
                    result_K.append(K)
                # if K_sp.Pic == 3 and K not in result_K_final:
                #     result_K.append(K.copy())
            elif counter == max_counter:
                if len(result_K) % 1000 == 0:
                    print(len(result_K))
                result_K.append((info_ridges, info_facets, queued_facets, counter, nbr_facets, nbr_ridges))
            elif nbr_facets <= 252:
                # We must find the set of every unclosed ridges and find the possible facets to close them for each of them
                facets_to_close_ridges = [None for ridge in ridges]
                unclosed_ridges = []
                for facet_index in queued_facets:
                    for ridge_data in G[facet_index][1]:
                        ridge_index = ridge_data[0]
                        if info_ridges[ridge_index] == 1:  # this means that the ridge is unclosed
                            unclosed_ridges.append(ridge_index)
                            # this is correct because the unclosed ridges are unique for every facets
                            # (if they were equal that would mean that the ridge is closed!)
                            facets_to_close_ridges[ridge_index] = []
                            for candidate_facet_index in ridge_data[1]:
                                if facet_layer[
                                    candidate_facet_index] == counter:  # this condition is enough, because the facets already used are in the previous layer
                                    facets_to_close_ridges[ridge_index].append(candidate_facet_index)
                            if facets_to_close_ridges[ridge_index] == []:
                                return graph_method_rec(info_ridges, info_facets, facet_layer, [], counter + 1,
                                                        nbr_facets, nbr_ridges, result_K,
                                                        max_counter)  # in this case, we found no candidate...

                # Now we have all the unclosed ridges with for each of them a list of candidates to close it
                # We need a recursive function to select a case
                unclosed_ridges.sort()  # like this, we can use dichotomie for finding an element in the list
                ridge_can_be_skipped = [False for unclosed_ridge in unclosed_ridges]
                list_next_cases = []
                current_case = []
                forbidden_facets = [False for facet in facets]

                def enumerate_cases_rec(k, l, current_case, ridge_can_be_skipped, forbidden_facets,
                                        transition_info_ridges, list_next_cases):
                    if k == len(unclosed_ridges):
                        if current_case not in list_next_cases:
                            list_next_cases.append(current_case)
                            new_info_ridges = transition_info_ridges
                            new_info_facets = info_facets.copy()
                            for facet_index in current_case:
                                new_info_facets[facet_index] = 1
                                # We update info ridges, the facets we add are supposed to close perfectly the last unclosed_ridges
                                # for ridge_data in G[facet_index][1]:
                                #     new_info_ridges[ridge_data[0]] += 1
                            graph_method_rec(new_info_ridges, new_info_facets, facet_layer, current_case, counter + 1,
                                             nbr_facets + len(current_case), nbr_ridges + len(unclosed_ridges),
                                             result_K, max_counter)
                    else:
                        if transition_info_ridges[unclosed_ridges[k]] == 2:
                            enumerate_cases_rec(k + 1, 0, current_case, ridge_can_be_skipped, forbidden_facets,
                                                transition_info_ridges, list_next_cases)
                        else:
                            if l + 1 < len(facets_to_close_ridges[unclosed_ridges[k]]):
                                enumerate_cases_rec(k, l + 1, current_case, ridge_can_be_skipped, forbidden_facets,
                                                    transition_info_ridges, list_next_cases)
                            facet_to_add = facets_to_close_ridges[unclosed_ridges[k]][l]
                            if not forbidden_facets[facet_to_add]:
                                new_ridge_can_be_skipped = ridge_can_be_skipped.copy()
                                new_forbidden_facets = forbidden_facets.copy()
                                new_transition_info_ridges = transition_info_ridges.copy()
                                for facet in facets_to_close_ridges[unclosed_ridges[k]]:
                                    new_forbidden_facets[facet] = True
                                problem = False
                                for ridge_data in G[facet_to_add][1]:
                                    ridge_index = ridge_data[0]
                                    new_transition_info_ridges[ridge_index] += 1
                                    if new_transition_info_ridges[ridge_index] > 2:
                                        problem = True
                                        break
                                    if info_ridges[ridge_index] == 1:
                                        position_in_unclosed_ridges = sc.dichotomie(unclosed_ridges, ridge_index)
                                        if position_in_unclosed_ridges > -1:
                                            new_ridge_can_be_skipped[position_in_unclosed_ridges] = True
                                            for facet in facets_to_close_ridges[
                                                unclosed_ridges[position_in_unclosed_ridges]]:
                                                new_forbidden_facets[facet] = True
                                if not problem:
                                    enumerate_cases_rec(k + 1, 0, current_case + [facet_to_add],
                                                        new_ridge_can_be_skipped, new_forbidden_facets,
                                                        new_transition_info_ridges,
                                                        list_next_cases)

                # print("hello")
                if [] not in facets_to_close_ridges:
                    enumerate_cases_rec(0, 0, current_case, ridge_can_be_skipped, forbidden_facets, info_ridges.copy(),
                                        list_next_cases)
            # print(len(list_next_cases))
            # if list_next_cases:
            #     for next_case in list_next_cases:

    facet_layer = [-1 for facet in facets]
    facet_layer[0] = 0
    find_layer(G, facet_layer, [0], 1)
    if starting_point == None:
        info_ridges = [0 for ridge in ridges]
        info_facets = [0 for facet in facets]
        queued_facets = [0]
        info_facets[0] = 1
        for data_ridges in G[0][1]:
            info_ridges[data_ridges[0]] += 1
        result_K_final = []
        graph_method_rec(info_ridges, info_facets, facet_layer, queued_facets, 1, 1, 11, result_K_final, max_counter)
        return result_K_final
    else:
        result_K_final = []
        (info_ridges, info_facets, queued_facets, counter, nbr_facets, nbr_ridges) = starting_point
        graph_method_rec(info_ridges, info_facets, facet_layer, queued_facets, counter, nbr_facets, nbr_ridges,
                         result_K_final, max_counter)
        return result_K_final


def initialize_graph_method3(facets, ridges, G, facet_layer, n, m):
    facets_for_ridges = []
    for ridge in ridges:
        facets_for_ridges.append([])
        for element in sc.list_2_pow[:m + 1]:
            if ridge | element != ridge:
                candidate_facet_for_ridge = ridge | element
                position = sc.dichotomie(facets, candidate_facet_for_ridge)
                if position != -1:
                    facets_for_ridges[-1].append(position)
    facets_for_ridges_with_layer = [
        [[position_facet for position_facet in facets_for_ridge if facet_layer[position_facet] == k] for
         facets_for_ridge in facets_for_ridges] for k in range(1, m - n + 1)]

    ridges_of_facets = [[ridge_data[0] for ridge_data in G[k][1]] for k in range(len(G))]
    print("graph method initialized")
    return facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets


def graph_method3_with_rec(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n, m,
                           starting_layer, max_layer,
                           starting_state=None):
    def graph_method_3_rec(info_facets, forbidden_facets, info_ridges, current_layer, results):
        if current_layer == max_layer:
            results.append((info_facets, forbidden_facets, info_ridges))
        # elif current_layer > m - n:
        #     if 1 not in info_ridges:
        #         results.append([facets[index] for index in range(len(info_facets)) if info_facets[index]])
        #     # K = sc.PureSimplicialComplex([facets[index] for index in range(len(info_facets)) if info_facets[index]])
        #     # if K.Pic == 4 and K.is_promising() and K.is_Z2_homology_sphere():
        #     #     results.append(K)
        #     # print([1*facet for facet in info_facets])
        else:
            def enumerate_cases(index_of_unclosed_ridges, k, l, current_info_facets, current_forbidden_facets,
                                current_info_ridges):
                if k == len(index_of_unclosed_ridges):
                    graph_method_3_rec(current_info_facets, current_forbidden_facets, current_info_ridges,
                                       current_layer + 1, result)
                else:
                    if current_info_ridges[index_of_unclosed_ridges[k]] == 2:
                        enumerate_cases(index_of_unclosed_ridges, k + 1, 0, current_info_facets,
                                        current_forbidden_facets,
                                        current_info_ridges)
                    else:  # if the ridge hasn't been closed yet
                        if facets_for_ridges[
                            index_of_unclosed_ridges[k]] != []:  # must be non empty
                            if l + 1 < len(facets_for_ridges[index_of_unclosed_ridges[k]]):
                                enumerate_cases(index_of_unclosed_ridges, k, l + 1, current_info_facets,
                                                current_forbidden_facets, current_info_ridges)
                            index_candidate_facet = \
                                facets_for_ridges[index_of_unclosed_ridges[k]][l]
                            if not forbidden_facets[index_candidate_facet]:  # if the facet is not forbidden
                                new_current_info_facets = current_info_facets.copy()
                                new_current_forbidden_facets = current_forbidden_facets.copy()
                                new_current_info_ridges = current_info_ridges.copy()
                                new_current_info_facets[index_candidate_facet] = True
                                for index_newly_closed_ridge in ridges_of_facets[index_candidate_facet]:
                                    new_current_info_ridges[index_newly_closed_ridge] += 1
                                    if new_current_info_ridges[index_newly_closed_ridge] == 2:
                                        for index_new_forbidden_facet in facets_for_ridges[index_newly_closed_ridge]:
                                            new_current_forbidden_facets[index_new_forbidden_facet] |= True
                                # new_current_info_ridges[ridges_of_facets[index_candidate_facet]] += 1
                                # for index_newly_closed_ridge in \
                                # np.where(np.logical_and(new_current_info_ridges == 2, current_info_ridges == 1))[0]:
                                #     new_current_forbidden_facets[
                                #         facets_for_ridges[int(index_newly_closed_ridge)]] = True  # changer
                                enumerate_cases(index_of_unclosed_ridges, k + 1, 0, new_current_info_facets,
                                                new_current_forbidden_facets, new_current_info_ridges)

            index_of_unclosed_ridges = [index_ridge for index_ridge in range(len(info_ridges)) if
                                        info_ridges[index_ridge] == 1]
            if index_of_unclosed_ridges != []:
                if current_layer<=m-n:
                    enumerate_cases(index_of_unclosed_ridges, 0, 0, info_facets.copy(), forbidden_facets.copy(),
                            info_ridges.copy())
            else:
                if max_layer == -1:
                    result.append([facets[index] for index in range(len(info_facets)) if info_facets[index]])
                    # print(result[-1])

    if starting_state == None:
        # here is some initialization
        info_facets = [False for facet in facets]
        forbidden_facets = [False for facet in facets]
        info_facets[0] = True
        forbidden_facets[0] = True
        info_ridges = [0 for ridge in ridges]
        for ridge_index in ridges_of_facets[0]:
            info_ridges[ridge_index] += 1
        result = []
        graph_method_3_rec(info_facets, forbidden_facets, info_ridges, starting_layer, result)
        return result
    else:
        info_facets, forbidden_facets, info_ridges = starting_state
        result = []
        graph_method_3_rec(info_facets, forbidden_facets, info_ridges, starting_layer, result)
        return (result)


def graph_method3_with_iter(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n, m,
                            starting_layer, max_layer,
                            starting_state=None):
    def graph_method_3_rec(info_facets, forbidden_facets, info_ridges, current_layer, results):
        if current_layer == max_layer:
            results.append((info_facets, forbidden_facets, info_ridges))
        # elif current_layer >= m - n:
        #     if 1 not in info_ridges:
        #         # results.append([facets[index] for index in range(len(info_facets)) if info_facets[index]])
        #         K = sc.PureSimplicialComplex([facets[index] for index in range(len(info_facets)) if info_facets[index]])
        #         if K.Pic == 2 and K.is_promising() and K.is_Z2_homology_sphere():
        #             results.append(K)
        #         print([1*facet for facet in info_facets])
        else:
            index_of_unclosed_ridges = [index_ridge for index_ridge in range(len(info_ridges)) if
                                        info_ridges[index_ridge] == 1]
            if index_of_unclosed_ridges != []:
                index_to_start = [0 for index_ridge in index_of_unclosed_ridges]
                info_facets_data = [info_facets.copy() for index_ridge in index_of_unclosed_ridges]
                info_ridges_data = [info_ridges.copy() for index_ridge in index_of_unclosed_ridges]
                forbidden_facets_data = [forbidden_facets.copy() for index_ridge in index_of_unclosed_ridges]
                k = 0
                l = 0
                going_forward = True
                while k >= 0:
                    # print(index_to_start)
                    if going_forward:
                        if k == len(index_of_unclosed_ridges):
                            graph_method_3_rec(info_facets_data[k - 1].copy(), forbidden_facets_data[k - 1].copy(),
                                               info_ridges_data[k - 1].copy(),
                                               current_layer + 1, result)
                            going_forward = False
                            k -= 1
                        else:
                            if k == 0:
                                info_ridges_data[k] = info_ridges.copy()
                            else:
                                info_ridges_data[k] = info_ridges_data[k - 1].copy()

                            if info_ridges_data[k][index_of_unclosed_ridges[k]] == 1:
                                if k == 0:
                                    forbidden_facets_data[k] = forbidden_facets.copy()
                                else:
                                    forbidden_facets_data[k] = forbidden_facets_data[k - 1].copy()

                                l = index_to_start[k]
                                if l < len(facets_for_ridges_with_layer[current_layer][index_of_unclosed_ridges[k]]):
                                    index_candidate_facet = \
                                        facets_for_ridges_with_layer[current_layer][index_of_unclosed_ridges[k]][l]
                                    if k == 0:
                                        info_facets_data[k] = info_facets.copy()
                                    else:
                                        info_facets_data[k] = info_facets_data[k - 1].copy()
                                    if not forbidden_facets_data[k][index_candidate_facet]:
                                        info_facets_data[k][index_candidate_facet] = True
                                        for index_newly_closed_ridge in ridges_of_facets[index_candidate_facet]:
                                            info_ridges_data[k][index_newly_closed_ridge] += 1
                                            if info_ridges_data[k][index_newly_closed_ridge] == 2:
                                                for index_new_forbidden_facet in facets_for_ridges[
                                                    index_newly_closed_ridge]:
                                                    forbidden_facets_data[k][index_new_forbidden_facet] |= True
                                    index_to_start[k] = l + 1
                                    k += 1
                                else:
                                    going_forward = False

                            else:
                                info_ridges_data[k] = info_ridges_data[k - 1].copy()
                                forbidden_facets_data[k] = forbidden_facets_data[k - 1].copy()
                                info_facets_data[k] = info_facets_data[k - 1].copy()
                                index_to_start[k] = l + 1
                                k += 1
                    else:
                        while (k > 0 and (
                                info_ridges_data[k - 1][index_of_unclosed_ridges[k]] == 2 or index_to_start[k] == len(
                            facets_for_ridges_with_layer[current_layer][index_of_unclosed_ridges[k]]))) or (
                                k == 0 and (info_ridges[index_of_unclosed_ridges[k]] == 2 or index_to_start[k] == len(
                            facets_for_ridges_with_layer[current_layer][index_of_unclosed_ridges[k]]))):
                            index_to_start[k] = 0

                            k -= 1
                        going_forward = True
            else:
                if max_layer == -1:
                    results.append([facets[index] for index in range(len(info_facets)) if info_facets[index]])
                    # K = sc.PureSimplicialComplex(
                    #     [facets[index] for index in range(len(info_facets)) if info_facets[index]])
                    # if K.Pic == 3 and K.is_promising() and K.is_Z2_homology_sphere():
                    #     results.append(K)

    if starting_state == None:
        # here is some initialization
        info_facets = [False for facet in facets]
        forbidden_facets = [False for facet in facets]
        info_facets[0] = True
        forbidden_facets[0] = True
        info_ridges = [0 for ridge in ridges]
        for ridge_index in ridges_of_facets[0]:
            info_ridges[ridge_index] += 1
        result = []
        graph_method_3_rec(info_facets, forbidden_facets, info_ridges, starting_layer, result)
        return result
    else:
        info_facets, forbidden_facets, info_ridges = starting_state
        result = []
        graph_method_3_rec(info_facets, forbidden_facets, info_ridges, starting_layer, result)
        return (result)


def construct_complete_graph(n, m):
    list_2_pow = [1]
    for k in range(m + 1):
        list_2_pow.append(list_2_pow[-1] * 2)
    full_simplex = list_2_pow[m] - 1
    facets = []
    for facet_iter in combinations(range(1, m + 1), n):
        facets.append(sc.face_to_binary(list(facet_iter), m))
    ridges = []
    for ridge_iter in combinations(range(1, m + 1), n - 1):
        ridges.append(sc.face_to_binary(list(ridge_iter), m))
    G = [(k, []) for k in range(len(facets))]
    for i in range(len(G)):
        for element1 in list_2_pow:
            if element1 | facets[i] == facets[i]:
                ridge = element1 ^ facets[i]
                position_ridge = sc.dichotomie(ridges, ridge)
                G[i][1].append((position_ridge, []))
                other_vertices = full_simplex ^ facets[i]
                for element2 in list_2_pow:
                    if element2 | other_vertices == other_vertices:
                        position_of_connected_facet = sc.dichotomie(facets, element2 | ridge)
                        if position_of_connected_facet >= 0:
                            G[i][1][-1][1].append(position_of_connected_facet)
    print("Graph constructed")
    return facets, ridges, G


def find_layer_v2(facets,ridges, m, starting_facet_index):
    list_2_pow = [1]
    for k in range(m + 1):
        list_2_pow.append(list_2_pow[-1] * 2)
    facets_layer = [-1 for facet in facets]
    current_layer = 0
    facets_layer[starting_facet_index] = current_layer
    current_layer += 1
    while -1 in facets_layer:
        for index_facet in range(len(facets)):
            if facets_layer[index_facet] == current_layer - 1:
                for index_neighbour in range(len(facets)):
                    if facets_layer[index_neighbour] == -1:
                        print(facets[index_neighbour] & facets[index_facet])
                        if facets[index_neighbour] & facets[index_facet] in ridges:
                            facets_layer[index_neighbour] = current_layer
        current_layer+=1
    return facets_layer
