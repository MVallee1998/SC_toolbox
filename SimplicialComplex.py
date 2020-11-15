import sqlite3
import numpy.polynomial as npp
import binary_search_tree as bst
import json
import datetime
from itertools import combinations, permutations
import timeit
import Z2_linear_algebra
import Betti_numbers as bnbr

list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1] * 2)


def face_to_binary(face, m):
    binary_face = m * [0]
    for k in face:
        binary_face[m - k] = 1
    return int(''.join(map(str, binary_face)), 2)


class PureSimplicialComplex:
    def __init__(self, facets=None, MNF_set=None, n=0):
        self.m = 0
        self.n = n
        self.facets = facets
        self.MNF_set = MNF_set
        if facets:
            self.facets = facets
            if type(facets[0]) == list:
                self.n = len(facets[0])
                labels = []
                for facet in facets:
                    for i in facet:
                        if i not in labels:
                            labels.append(i)
                self.m = len(labels)
                self.list_2_pow = [1]
                for k in range(16):
                    self.list_2_pow.append(self.list_2_pow[-1] * 2)
                self.facets_bin = [face_to_binary(facet, self.m) for facet in self.facets]
            else:
                self.n = 0
                self.m = 0
                self.facets_bin = facets
                self.facets_bin.sort()
                last_facet = self.facets_bin[-1]
                two_pow = 1
                while last_facet % two_pow != last_facet:
                    if two_pow | last_facet == last_facet:
                        self.n += 1
                    self.m += 1
                    two_pow *= 2
                self.list_2_pow = [1]
                for k in range(16):
                    self.list_2_pow.append(self.list_2_pow[-1] * 2)
        elif MNF_set:
            if not n:
                raise Exception
            else:
                self.MNF_set = MNF_set
                self.n = n
                labels = []
                for MNF in MNF_set:
                    for i in MNF:
                        if i not in labels:
                            labels.append(i)
                self.m = len(labels)
                self.MNF_set_bin = [face_to_binary(MNF, self.m) for MNF in self.MNF_set]
                self.compute_facets_from_MNF_set()
        self.FP_bin = None
        self.Pic = self.m - self.n
        self.f_vector = None
        self.h_vector = None
        self.g_vector = None

        self.H = None
        self.unclosed_ridges = None
        self.closed_faces = None

    def create_FP(self):
        if not self.FP_bin:
            self.FP_bin = [[] for i in range(self.n)]
            faces_set = [bst.Node() for k in range(self.n)]
            facet_data = [(facet, []) for facet in self.facets_bin]
            faces_set[-1].insertList(facet_data)
            for k in range(self.n - 2, -1, -1):  # dimension
                faces = []
                faces_set[k + 1].TreeToList(faces)
                self.FP_bin[k + 1] = faces.copy()
                for l in range(len(faces)):  # treat every face of dimension k
                    face = faces[l][0]
                    for element in self.list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                        if element | face == face:
                            subface = element ^ face
                            faces_set[k].insert((subface, [l]))
            vertices = []
            faces_set[0].TreeToList(vertices)
            self.FP_bin[0] = vertices.copy()

    def create_f_vector(self):
        if not self.FP_bin:
            self.create_FP()
        self.f_vector = [1]
        for i in range(self.n):
            self.f_vector.append(self.FP_bin[i].size())

    def create_h_vector(self):
        if not self.f_vector:
            self.create_f_vector()
        f_vector_copy = list(self.f_vector)
        f_vector_copy.reverse()
        F_pol = npp.Polynomial(f_vector_copy)
        P = npp.Polynomial([-1, 1])
        H_pol = F_pol(P)
        self.h_vector = list(H_pol)
        self.h_vector = [int(h_i) for h_i in self.h_vector]

    def create_g_vector(self):
        if not self.h_vector:
            self.create_h_vector()
        self.g_vector = []
        for i in range(1, int(self.n / 2) + 1):
            self.g_vector.append(self.h_vector[i] - self.h_vector[i - 1])

    def faces_set_to_MNF_set(self):
        temp_faces_list = [self.list_2_pow[self.m] - 1]
        k = self.m
        subface_set = bst.Node()
        subface_set.insertList(temp_faces_list)
        # we enumerate all the faces until the time we will have the right dimension faces
        while k > self.n:
            subface_set = bst.Node()
            for face in temp_faces_list:  # treat every face of dimension k
                for element in self.list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face == face:
                        subface = element ^ face
                        if not subface_set.findval(subface):
                            if k > self.n or (k == self.n and (not dichotomie(subface, self.FP_bin[k - 1]))):
                                subface_set.insert((subface, []))
            temp_faces_list = []
            subface_set.TreeToList(temp_faces_list)
            k = k - 1
        # we then create the set with all subsets of [m] of size <= n not in  the given faces set
        all_non_faces = [bst.Node() for k in range(self.n + 1)]
        all_non_faces[-1] = subface_set
        for k in range(self.n - 1, -1, -1):  # dimension
            faces = []
            all_non_faces[k + 1].TreeToList(faces)
            for l in range(len(faces)):  # treat every face of dimension k
                face = faces[l]
                for element in self.list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face == face:
                        subface = element ^ face
                        if (not (all_non_faces[k].findval(subface))) and (not dichotomie(subface, self.FP_bin[k])):
                            all_non_faces[k].insert((subface, [l + 1]))
        MNF_set = bst.Node()
        # Here we will try every face of dimension k
        for k in range(1, self.n + 1):
            faces_to_test = []
            all_non_faces[k].TreeToList(faces_to_test)
            for face_to_test in faces_to_test:
                is_MNF = True
                for element in self.list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face_to_test == face_to_test:
                        subface = element ^ face_to_test
                        if not (dichotomie(subface, self.FP_bin[k - 1])):
                            is_MNF = False
                            break
                if is_MNF:
                    MNF_set.insert((face_to_test, [1]))
        self.MNF_set_bin = []
        MNF_set.TreeToList(self.MNF_set_bin)
        self.MNF_set_bin = self.MNF_set_bin[:][0]

    def compute_facets_from_MNF_set(self):  # TO CHANGE
        M = range(1, self.m + 1)
        candidate_facets_iter = combinations(M, self.n)
        candidate_facets = []
        for facet_iter in candidate_facets_iter:
            candidate_facets.append(sum([self.list_2_pow[k] for k in facet_iter]))
        self.facets_bin = []
        for facet in candidate_facets:
            is_a_facet = True
            for MNF in self.MNF_set:  # We check if it the facet does not contain some minimal non-face
                if MNF | facet == facet:
                    is__a_facet = False
                    break
            if is_a_facet:
                self.facets_bin.append(facet)

    def binary_to_face(self, face_bin):
        face = []
        for k in range(self.m):
            if face_bin | self.list_2_pow[k] == face_bin:
                face.append(k)
        return face

    def is_closed(self, F=None):
        if F:
            Link_of_F = Link_of(self, F)
            ridges_set = bst.Node()
            if Link_of_F.n < 2: return True
            for k in range(len(Link_of_F.facets_bin)):
                facet = Link_of_F.facets_bin[k]
                for element in Link_of_F.list_2_pow:
                    if element | facet == facet:
                        ridges_set.insert((element ^ facet, [k]))
            closed = True
            ridges_data_list = []
            ridges_set.TreeToList(ridges_data_list)
            for ridge_data in ridges_data_list:
                if len(ridge_data[1]) != 2:
                    closed = False
                    break
            return closed
        else:
            if not self.FP_bin:
                ridges_set = bst.Node()
                if self.n < 2: return True
                for k in range(len(self.facets_bin)):
                    facet = self.facets_bin[k]
                    for element in self.list_2_pow:
                        if element | facet == facet:
                            ridges_set.insert((element ^ facet, [k]))
                ridges_data_list = []
                ridges_set.TreeToList(ridges_data_list)
            else:
                ridges_data_list = self.FP_bin[self.n - 2]
            closed = True
            for ridge_data in ridges_data_list:
                if len(ridge_data[1]) != 2:
                    closed = False
                    break
            return closed

    def list_closed_faces(self):
        if not self.FP_bin:
            self.create_FP()
        if self.closed_faces == None:
            self.closed_faces = []
            for k in range(self.n - 2):
                for face_data in self.FP_bin[k]:
                    if self.is_closed(face_data[0]):
                        self.closed_faces.append(face_data[0])

    def list_unclosed_ridges(self):
        if not self.unclosed_ridges:
            ridges_set = bst.Node()
            if self.n < 2: return []
            for k in range(len(self.facets_bin)):
                facet = self.facets_bin[k]
                for element in self.list_2_pow:
                    if element | facet == facet:
                        ridges_set.insert((element ^ facet, [k + 1]))
            ridges_data_list = []
            ridges_set.TreeToList(ridges_data_list)
            unclosed_ridges = []
            for ridge_data in ridges_data_list:
                if len(ridge_data[1]) != 2:
                    unclosed_ridges.append(ridge_data[0])
            unclosed_ridges.sort()
            self.unclosed_ridges = unclosed_ridges

    def is_promising(self):
        if not self.FP_bin:
            self.create_FP()
        self.list_closed_faces()
        for closed_face in self.closed_faces:
            Link_K_of_F = Link_of(self, closed_face)
            Link_K_of_F.Z2_Betti_numbers()
            if not Link_K_of_F.is_Z2_homology_sphere():
                return False
        return True

    def is_candidate(self, S, aimed_m):
        if S <= self.facets_bin[-1]:
            return False
        boundary_of_S = []
        self.list_closed_faces()
        for element in self.list_2_pow:
            if element | S == S:
                boundary_of_S.append(element ^ S)
        boundary_of_S.sort()
        self.list_unclosed_ridges()
        unclosed_ridges = self.unclosed_ridges
        for ridge in boundary_of_S:
            if dichotomie(unclosed_ridges, ridge) >= 0:
                for closed_face in self.closed_faces:
                    if S | closed_face == S:
                        return False
                return True
        return False

    def Z2_Betti_numbers(self):
        if not self.H:
            if not self.FP_bin:
                self.create_FP()
            FP = [[face_data[0] for face_data in self.FP_bin[k]] for k in range(self.n)]
            boundary_matrices = [[] for k in range(self.n)]
            boundary_matrices[0] = [[] for k in range(self.m)]
            for k in range(1, self.n):
                for face in FP[k]:
                    boundary_matrices[k].append([])
                    for element in self.list_2_pow:
                        if face | element == face:
                            boundary_matrices[k][-1].append(dichotomie(FP[k - 1], face ^ element))
                    boundary_matrices[k][-1].sort()
            # im = 0
            # H = []
            # for i in range(1, self.n):
            #     boundary_matrix = Z2_linear_algebra.Z2Array(len(boundary_matrices[i]), boundary_matrices[i - 1])
            #     ker = boundary_matrix.n - boundary_matrix.Z2_rank()-1
            #     H.append(ker - im)
            #     im = boundary_matrix.n - ker
            # H.append(len(boundary_matrices[self.n - 1]) - im)
            self.H = bnbr.computeBettiNumbers(boundary_matrices)

    def is_Z2_homology_sphere(self):
        if not self.H:
            self.Z2_Betti_numbers()
        if self.H == [2]:
            return True
        L = [1]
        for k in range(1, self.n - 1):
            L.append(0)
        L.append(1)
        return self.H == L

    def is_minimal_lexico_order(self):
        closed_vertices = []
        for v in range(self.m):
            link_of_v = Link_of(self, v)
            if self.is_closed(self.list_2_pow[v]):
                closed_vertices.append(v)

        ridges = []
        minimal_facets_bin = self.facets_bin.copy()
        for v in closed_vertices:
            ridges_containing_v = [ridge_data[0] for ridge_data in self.FP_bin[self.n - 2] if
                                   ridge_data[0] | self.list_2_pow[v] == ridge_data]
            for ridge in ridges_containing_v:
                ridge_labels_to_modify = self.binary_to_face(ridge ^ self.list_2_pow[v])
                for ridge_labels_permu_iter in permutations(ridge_labels_to_modify):
                    F = []
                    for facet in self.facets_bin:
                        if facet | ridge == facet:
                            F += (self.binary_to_face(facet ^ ridge))
                    remaining_labels = self.binary_to_face((self.list_2_pow[self.m] - 1) ^ ridge ^ F[0] ^ F[1])
                    for F_perm_iter in permutations(F):
                        for remaining_labels_perm_iter in permutations(remaining_labels):
                            old_labels = [v]
                            old_labels += list(ridge_labels_permu_iter)
                            old_labels += list(F_perm_iter)
                            old_labels += list(remaining_labels_perm_iter)
                            relabeled_facets = relabel_facets(self, old_labels)
                            if relabeled_facets < minimal_facets_bin:
                                return False
        return True


def relabel_facets(K, old_labels):
    if len(old_labels) != K.m:
        raise Exception
    new_facets = [0] * len(K.facets_bin)
    for k in range(len(new_facets)):
        for l in range(K.m):
            if K.facets_bin[k] ^ K.list_2_pow[old_labels[l]] == K.facets_bin[k]:
                # the case when 01 or 10
                new_facets[k] += K.list_2_pow[l]
    new_facets.sort()
    return new_facets


def Link_of(K, F):
    if not K.FP_bin:
        K.create_FP()
    k = 0
    l = 0
    while F % K.list_2_pow[l] != F:
        if K.list_2_pow[l] | F == F:
            k += 1
        l += 1
    is_a_face = False
    for facet in K.facets_bin:
        if F | facet == facet:
            is_a_face = True
            break
    if not is_a_face:
        raise Exception
    facets_of_Link = []
    complementary_faces = [facet_data[0] for facet_data in K.FP_bin[K.n - k - 1]]
    for complementary_face in complementary_faces:
        if dichotomie(K.facets_bin, complementary_face ^ F) >= 0:
            facets_of_Link.append(complementary_face)
    unused_labels = []
    for i in range(K.m):
        unused = True
        for facet in facets_of_Link:
            if K.list_2_pow[i] | facet == facet:
                unused = False
                break
        if unused:
            unused_labels.append(i)
    for i in range(len(unused_labels) - 1, -1, -1):
        for j in range(len(facets_of_Link)):
            unshifted_bits = facets_of_Link[j] % K.list_2_pow[unused_labels[i]]
            facets_of_Link[j] = ((facets_of_Link[j] ^ unshifted_bits) >> 1) ^ unshifted_bits
    if len(unused_labels) == 0:
        return PureSimplicialComplex(facets_of_Link)
    m_of_link = K.m
    label_unused = True
    while m_of_link > 0 and label_unused:
        for facet in facets_of_Link:
            if facet ^ K.list_2_pow[m_of_link - 1] == facet:
                label_unused = False
                break
        m_of_link -= 1
    return PureSimplicialComplex(facets_of_Link)


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip().decode('ASCII') for x in data]
    return data


def give_dim(facets):
    n = len(facets[0])
    m = max([facet[-1] for facet in facets])
    pic = m - n
    return n, m, pic


def dicho_search(x, l, a, b):
    if b - a <= 1:
        return l[a] == x or l[b] == x
    else:
        i = int((a + b) / 2)
        if x > l[i]:
            return dicho_search(x, l, i, b)
        elif x < l[i]:
            return dicho_search(x, l, a, i)
        else:
            return True


def dichotomie(t, v):
    a = 0
    b = len(t) - 1
    while a <= b:
        m = (a + b) // 2
        if t[m] == v:
            # on a trouvé v
            return m
        elif t[m] < v:
            a = m + 1
        else:
            b = m - 1
    # on a a > b
    return -1


def Garrison_Scott(K):
    facets = K.facets
    K.create_FP()
    list_char_funct = []
    '''we create all the non zero elements of Z_2^n'''
    list_2_pow = [1]
    for k in range(1, K.m + 1):
        list_2_pow.append(list_2_pow[-1] * 2)
    list_S = []
    for k in range(K.n, K.m):
        list_S.append(list(range(1, K.list_2_pow[K.n])))
    K_full_FP = []
    for k in range(K.n):
        for face in K.FP_bin[k]:
            K_full_FP.append(face)
    current_lambda = []
    i = K.n
    for k in range(K.n):
        current_lambda.append(list_2_pow[k])
    for k in range(K.n, K.m):
        current_lambda.append(0)
    while i >= K.n:
        if len(list_char_funct) % 1000 == 0:
            print(len(list_char_funct))
        for face_bin in K_full_FP:
            if (face_bin < list_2_pow[i + 1]) and ((face_bin | list_2_pow[i]) == face_bin):
                linear_comb = 0
                face_bin_transformed = face_bin ^ list_2_pow[i]
                for k in range(K.m):
                    if face_bin_transformed | list_2_pow[k] == face_bin_transformed:
                        linear_comb = current_lambda[k] ^ linear_comb
                index = dichotomie(list_S[K.n - i], linear_comb)
                if index != -1:
                    list_S[K.n - i].pop(index)
        while len(list_S[K.n - i]) > 0:
            current_lambda[i] = list_S[K.n - i].pop()
            if i + 1 == K.m:
                list_char_funct.append(current_lambda.copy())
            else:
                i += 1
                break
        if len(list_S[K.n - i]) == 0:
            list_S[K.n - i] = list(range(1, list_2_pow[K.n]))
            i -= 1
    return list_char_funct


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
        candidate_facets.append((list_2_pow[m] - 1) ^ face_to_binary(cofacet, m))
    candidate_facets.sort()
    ridges = []
    for facet in candidate_facets:
        for element in list_2_pow:
            if element | facet == facet:
                ridge = element ^ facet
                if not ridge in ridges:
                    ridges.append(ridge)
    ridges.sort()
    return candidate_facets, ridges


def construct_graph(char_funct, n, m):
    full_simplex = list_2_pow[m] - 1
    facets, ridges = enumerate_facets_and_ridges(char_funct, n, m)
    G = [(k, []) for k in range(len(facets))]
    for i in range(len(G)):
        for element1 in list_2_pow:
            if element1 | facets[i] == facets[i]:
                ridge = element1 ^ facets[i]
                position_ridge = dichotomie(ridges, ridge)
                G[i][1].append((position_ridge, []))
                other_vertices = full_simplex ^ facets[i]
                for element2 in list_2_pow:
                    if element2 | other_vertices == other_vertices:
                        position_of_connected_facet = dichotomie(facets, element2 ^ ridge)
                        if position_of_connected_facet >= 0:
                            G[i][1][-1][1].append(position_of_connected_facet)

        # for j in range(len(facets)):
        #     ridge_index = dichotomie(ridges, facets[i] & facets[j])
        #     if ridge_index >= 0:
        #         G[i][2].append((j,ridge_index))
        #         if ridge_index not in G[i][1]:
        #             G[i][1].append(ridge_index)
    return facets, ridges, G


def enumerate_cases(list_of_max, k, current_list, results):
    if k == len(list_of_max):
        results.append(current_list[:len(list_of_max)])
    else:
        if list_of_max[k] == -1:
            enumerate_cases(list_of_max, k + 1, current_list + [-1], results)
        else:
            enumerate_cases(list_of_max, k + 1, current_list + [0], results)
            if current_list[k] < list_of_max[k] - 1:
                new_list = current_list.copy()
                new_list[k] += 1
                enumerate_cases(list_of_max, k, new_list, results)


def graph_method(facets, ridges, G):
    def graph_method_rec(info_ridges_arg, info_facets_arg, queued_facets, counter):
        info_ridges = info_ridges_arg.copy()
        info_facets = info_facets_arg.copy()
        if counter >= 2:
            K = []
            for i in range(len(info_facets)):
                if info_facets == 1:
                    K.append(facets[i])
        if not queued_facets:
            if 1 not in info_ridges:
                K = []
                for i in range(len(info_facets)):
                    if info_facets == 1:
                        K.append(facets[i])
                print("coucou",K)
        else:
            list_cases = []
            for facet_index in queued_facets:
                list_of_max = []
                for list_neighbours in G[facet_index][1]:
                    position_ridge = list_neighbours[0]
                    if info_ridges[position_ridge] == 2:
                        list_of_max.append(-1)
                    else:
                        list_of_max.append(len(list_neighbours[1]))
                list_cases_current = []
                enumerate_cases(list_of_max, 0, [0], list_cases_current)
                list_cases.append(list_cases_current.copy())
            list_cases_ref = list_cases.copy()
            info_ridges_ref = info_ridges.copy()
            info_facets_ref = info_facets.copy()
            # print(sum([len(list_cases[i])for i in range(len(list_cases))]))
            k = 0
            list_facets_to_add = []
            while k>-1:
                problem_ridges = False
                if list_cases[k] == []:
                    list_cases[k] = list_cases_ref[k].copy()
                    k -= 1
                    continue
                if k == 0:
                    list_facets_to_add = []
                    info_ridges = info_ridges_ref.copy()
                    info_facets = info_facets_ref.copy()
                index_facets_to_add = list_cases[k].pop()
                for l in range(len(index_facets_to_add)):
                    index_facet = index_facets_to_add[l]
                    if info_facets[G[queued_facets[k]][1][l][1][index_facet]] != 1:
                        list_facets_to_add.append(G[queued_facets[k]][1][l][1][index_facet])
                        info_facets[G[queued_facets[k]][1][l][1][index_facet]] = 1
                        for ridge_data in G[G[queued_facets[k]][1][l][1][index_facet]][1]:
                            info_ridges[ridge_data[0]] += 1
                            if info_ridges[ridge_data[0]] > 2: # maybe this problem
                                problem_ridges = True
                                break
                        if problem_ridges:
                            break
                if problem_ridges:
                    break
                if (k + 1 == len(list_cases)):
                    graph_method_rec(info_ridges, info_facets, list_facets_to_add,counter+1)
                else:
                    k+=1
    info_ridges = [0 for ridge in ridges]
    info_facets = [0 for facet in facets]
    queued_facets = [0]
    info_facets[0] = 1
    for data_ridges in G[0][1]:
        info_ridges[data_ridges[0]] = 1

    graph_method_rec(info_ridges, info_facets, queued_facets, 0)


#
def Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m):
    if pile == []:
        return results
    K = pile.pop()
    print(K.facets_bin, results)
    if not K.is_closed():
        print("coucou")
        for face in candidate_facets_ref:
            if K.is_candidate(face, aimed_m):
                new_K = PureSimplicialComplex(K.facets_bin + [face])
                pile.append(new_K)
        return Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)
    else:
        results.append(K)
        return Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)


def Choi_Vallee_algo(pile, candidate_facets, info_facets, ridges, info_ridges):
    return True


facets, ridges, G = construct_graph([3,5,6,7,4, 2, 1], 4,7)
graph_method(facets, ridges, G)

# def find_minimal_facets_set(facets_set):
#     facets_set.sort()
#     n = len(facets_set[0])
#     # if the label 1 is not used, we must change it to 1 and change every other labels
#     if facets_set[0][0] != 1:
#         a = facets_set[0][0]
#         for i in range(len(facets_set)):
#             for j in range(n):
#                 facets_set[i][j] -= (a - 1)
#     m = max([facet[-1] for facet in facets_set])
#     pic = m - n
#     first_supposed_facet = list(range(1, n + 1))
#     # if the first facet is not the regular n-simplex [1,...,n] then we relabel every faces to match it
#     if facets_set[0] != first_supposed_facet:
#         permutation = facets_set[0]
#         for i in range(len(facets_set)):
#             for j in range(n):
#                 a = facets_set[i][j]
#                 if a <= n:
#                     facets_set[i][j] = permutation[a]
#             facets_set[i].sort()
#         facets_set.sort()
#     return (facets_set)
#     # At this point, only the labels (n+1, ..., m) should be modified to obtain the minimal pure simplicial complex
#     # for the lexicographic order
#     minimal_facets_set = [facet.copy() for facet in facets_set]
#     for last_labels_permutation_iter in permutations(range(n + 1, m + 1)):
#         temp_facets_set = [facet.copy() for facet in facets_set]
#         last_labels_permutation = list(last_labels_permutation_iter)
#         for i in range(len(facets_set)):
#             for j in range(n):
#                 a = temp_facets_set[i][j]
#                 if a > n:
#                     temp_facets_set[i][j] = last_labels_permutation[a - n - 1]
#             temp_facets_set[i].sort()
#         temp_facets_set.sort()
#         if temp_facets_set < minimal_facets_set:
#             minimal_facets_set = [facet.copy() for facet in temp_facets_set]
#     return minimal_facets_set


# facets_list = read_file('./PLS_6_3')
#
# start = timeit.default_timer()
# for k in range(len(facets_list)):
#     facets = facets_list[k]
#     list_lambdas = Garrison_Scott(json.loads(facets))
#     if len(list_lambdas) == 1:
#         print(facets)
# stop = timeit.default_timer()
# print(stop - start)
