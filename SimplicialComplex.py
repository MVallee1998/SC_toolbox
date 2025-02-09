import json
from tqdm import tqdm
from itertools import combinations, permutations

import numpy as np
from functools import reduce
import numpy.polynomial as npp

import Betti_numbers as bnbr
import SimplicialComplex
import Z2_linear_algebra as Z2la

# sys.setrecursionlimit(1000)

G_vector = [2, 6, 10, 20, 30, 50, 70, 105, 140, 196, 252]
list_2_pow = [1]
for k in range(30):
    list_2_pow.append(list_2_pow[-1] * 2)


def face_to_binary(face, m):
    binary_face = m * [0]
    for k in face:
        binary_face[m - k] = 1
    return int(''.join(map(str, binary_face)), 2)


def binary_to_face(x, m):
    return [k + 1 for k in range(m) if (list_2_pow[k] | x) == x]

def binary_to_face_0(x, m):
    return [k for k in range(m) if (list_2_pow[k] | x) == x]


class PureSimplicialComplex:
    def __init__(self, facets=None, MNF_set=None, n=0, FP_bin=None):
        self.m = 0
        self.n = n
        self.facets = facets
        self.MNF_set = MNF_set
        self.NF_set_bin = None
        self.facets_bin = []
        self.FP_bin = FP_bin
        self.MNF_set_bin = None
        self.f_vector = None
        self.h_vector = None
        self.g_vector = None
        self.H = None
        self.unclosed_ridges = None
        self.closed_faces = None
        self.list_faces = None
        if type(facets) == np.ndarray or type(facets)==list:
            self.facets = facets
            if type(facets[0]) == list:
                self.n = max([len(facet) for facet in self.facets])
                labels = []
                for facet in facets:
                    for i in facet:
                        if i not in labels:
                            labels.append(i)
                self.m = len(labels)
                self.facets_bin = np.array([face_to_binary(facet, self.m) for facet in self.facets])
            else:
                self.n = max([bin(MF).count("1") for MF in facets])
                self.m = bin(np.bitwise_or.reduce(np.array(facets))).count("1")
                np_facets = np.array(facets)
                or_facets = np.bitwise_or.reduce(np_facets)
                l = 0
                while l != self.m:
                    l = 0
                    while list_2_pow[l] | or_facets == or_facets:
                        l += 1
                    np_facets = (np_facets >> (l + 1)) << (l) | (((np_facets >> (l + 1)) << (l + 1)) ^ np_facets)
                    or_facets = np.bitwise_or.reduce(np_facets)
                # print(or_facets)
                self.facets_bin = np_facets.tolist()
                self.facets_bin.sort()
                self.facets = []
                for facet_bin in self.facets_bin:
                    self.facets.append(binary_to_face(facet_bin, self.m))
                self.facets.sort()
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
                self.facets = []
                for facet_bin in self.facets_bin:
                    self.facets.append(binary_to_face(facet_bin, self.m))
                self.facets.sort()
        self.Pic = self.m - self.n

    def create_FP(self):
        if not type(self.facets_bin):
            self.compute_facets_from_MNF_set()
        if not self.FP_bin:
            self.FP_bin = [dict() for i in range(self.n+1)]
            for facet_bin in self.facets_bin:
                self.FP_bin[bin(facet_bin).count("1") - 1][facet_bin] = []
            # self.FP_bin[-1] = dict.fromkeys(self.facets_bin, [])
            for k in range(self.n - 1, -1, -1):
                for face in self.FP_bin[k + 1]:
                    for element in list_2_pow[:self.m]:
                        if element | face == face:
                            subface = element ^ face
                            if subface not in self.FP_bin[k]:
                                self.FP_bin[k][subface] = [face]
                            else:
                                self.FP_bin[k][subface].append(face)

    def create_f_vector(self):
        if not self.FP_bin:
            self.create_FP()
        self.f_vector = [1]
        for i in range(self.n):
            self.f_vector.append(len(self.FP_bin[i]))

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

    def compute_NF_set(self):
        if not self.NF_set_bin:
            self.create_FP()
            # we enumerate all the faces until the time we will have the right dimension faces
            all_non_faces = [dict() for k in range(self.n + 1)]
            all_non_faces[-1] = dict.fromkeys(
                [face_to_binary(list(facet_iter), self.m) for facet_iter in
                 combinations(range(1, self.m + 1), self.n + 1)])
            for k in range(self.n-1 , -1, -1):  # dimension
                for face in all_non_faces[k + 1]:  # treat every face of dimension k
                    for element in list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                        if element | face == face:
                            subface = element ^ face
                            if subface not in self.FP_bin[k] and subface not in all_non_faces[k]:
                                all_non_faces[k][subface] = True
            self.NF_set_bin = all_non_faces.copy()

    def compute_MNF_set(self):
        self.compute_NF_set()
        self.MNF_set_bin = []
        # Here we will try every face of dimension k
        for k in range(1, self.n + 1):
            # print([binary_to_face(non_face,self.m) for non_face in all_non_faces[k]])
            for non_face_to_test in self.NF_set_bin[k]:
                is_MNF = True
                for element in list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | non_face_to_test == non_face_to_test:
                        subface = element ^ non_face_to_test
                        if subface not in self.FP_bin[k - 1]:
                            is_MNF = False
                            break
                if is_MNF and non_face_to_test not in self.MNF_set_bin:
                    self.MNF_set_bin.append(non_face_to_test)
        self.MNF_set_bin.sort()

    def MNF_bin_to_MNF(self):
        if not self.MNF_set:
            self.MNF_set = []
            for MNF_bin in self.MNF_set_bin:
                self.MNF_set.append(binary_to_face(MNF_bin, self.m))
            self.MNF_set.sort()

    def compute_facets_from_MNF_set(self):
        M = range(1, self.m + 1)
        candidate_facets_iter = combinations(M, self.n)
        candidate_facets = []
        for facet_iter in candidate_facets_iter:
            candidate_facets.append(sum([list_2_pow[k - 1] for k in facet_iter]))
        self.facets_bin = []
        print(candidate_facets)
        for facet in candidate_facets:
            is_a_facet = True
            for MNF in self.MNF_set_bin:  # We check if it the facet does not contain some minimal non-face
                if MNF | facet == facet:
                    is_a_facet = False
                    break
            if is_a_facet:
                self.facets_bin.append(facet)
        self.facets_bin.sort()
        print(self.facets_bin)
        print('coucou')
        # if self.facets_bin[0] != list_2_pow[self.n]-1:
        #     list_others = []
        #     missed = 0
        #     for v in range(self.m):
        #         if list_2_pow[v] | self.facets_bin[0] != self.facets_bin[0]:
        #     for v in range(self.n):
        #         if list_2_pow[v] | self.facets_bin[0] != self.facets_bin[0]:
        #             print("hello",v)
        #             for l in range(len(self.facets_bin)):
        #                 if (self.facets_bin[l]|list_2_pow[v]==self.facets_bin[l]) ^ (self.facets_bin[l]|list_2_pow[self.n+marker]==self.facets_bin[l]):
        #                     print("hello1")
        #                     modifier = (list_2_pow[v] + list_2_pow[self.n+marker])
        #                     self.facets_bin[l] = modifier ^ self.facets_bin[l]
        #             marker+=1
        #     self.facets_bin.sort()
        #     print(self.facets_bin)

    def binary_to_face(self, face_bin):
        face = []
        for k in range(self.m):
            if face_bin | list_2_pow[k] == face_bin:
                face.append(k)
        return face

    def is_closed(self, F=None):
        if F:
            Link_of_F = Link_of(self, F)
            ridges_set = dict()
            if Link_of_F.n < 2: return True
            for facet in Link_of_F.facets_bin:
                for element in list_2_pow:
                    if element | facet == facet:
                        subface = element ^ facet
                        if subface not in ridges_set:
                            ridges_set[subface] = [facet]
                        else:
                            ridges_set[subface].append(facet)
            for ridge_parents in ridges_set.values():
                if len(ridge_parents) != 2:
                    return False
            return True
        else:
            if not self.FP_bin:
                ridges_set = dict()
                if self.n < 2: return True
                for facet in self.facets_bin:
                    for element in list_2_pow:
                        if element | facet == facet:
                            subface = element ^ facet
                            if subface not in ridges_set:
                                ridges_set[subface] = [facet]
                            else:
                                ridges_set[subface].append(facet)
            else:
                ridges_set = self.FP_bin[self.n - 2]
            for ridge_data in ridges_set.items():
                if len(ridge_data[1]) != 2:
                    return False
            return True

    def list_closed_faces(self):
        self.create_FP()
        if self.closed_faces == None:
            self.closed_faces = []
            for k in range(self.n - 2):
                for face in self.FP_bin[k]:
                    if self.is_closed(face):
                        self.closed_faces.append(face)
            self.closed_faces.sort()

    def list_unclosed_ridges(self):
        if not self.unclosed_ridges:
            ridges_set = dict()
            if self.n < 2: return True
            for facet in self.facets_bin:
                for element in list_2_pow:
                    if element | facet == facet:
                        subface = element ^ facet
                        if subface not in ridges_set:
                            ridges_set[subface] = [facet]
                        else:
                            ridges_set[subface].append(facet)
            unclosed_ridges = []
            for ridge_data in ridges_set.items():
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
        for element in list_2_pow:
            if element | S == S:
                boundary_of_S.append(element ^ S)
        boundary_of_S.sort()
        self.list_unclosed_ridges()
        unclosed_ridges = self.unclosed_ridges
        unclosed_ridges.sort()
        for ridge in boundary_of_S:
            if ridge in unclosed_ridges:
                for closed_face in self.closed_faces:
                    if S | closed_face == S:
                        return False
                return True
        return False

    def Z2_Betti_numbers(self):
        if not self.H:
            self.create_FP()
            FP = [sorted([face for face in self.FP_bin[k]]) for k in range(self.n)]
            if len(self.facets_bin)==0:
                return []
            boundary_matrices = [[] for k in range(self.n)]
            boundary_matrices[0] = [[] for k in range(self.m)]
            for k in range(1, self.n):
                for face in FP[k]:
                    boundary_matrices[k].append([])
                    for element in list_2_pow:
                        if face | element == face:
                            boundary_matrices[k][-1].append(dichotomie(FP[k - 1], face ^ element))
                    boundary_matrices[k][-1].sort()
            self.H = bnbr.computeBettiNumbers(boundary_matrices)

    def is_Z2_homology_sphere(self):
        self.Z2_Betti_numbers()
        if self.n == 1 and self.H == [2]:
            return True
        L = [1]
        for k in range(1, self.n - 1):
            L.append(0)
        L.append(1)
        return self.H == L

    def is_minimal_lexico_order(self):
        closed_vertices = []
        for v in range(self.m):
            if self.is_closed(list_2_pow[v]):
                closed_vertices.append(v)
        minimal_facets_bin = self.facets_bin.copy()
        for v in closed_vertices:
            ridges_containing_v = [ridge for ridge in self.FP_bin[self.n - 2] if
                                   ridge | list_2_pow[v] == ridge]
            for ridge in ridges_containing_v:
                ridge_labels_to_modify = self.binary_to_face(ridge ^ list_2_pow[v])
                for ridge_labels_permu_iter in permutations(ridge_labels_to_modify):
                    F = []
                    for facet in self.facets_bin:
                        if facet | ridge == facet:
                            F += (self.binary_to_face(facet ^ ridge))
                    remaining_labels = self.binary_to_face(
                        (list_2_pow[self.m] - 1) ^ ridge ^ list_2_pow[F[0]] ^ list_2_pow[F[1]])
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

    def find_minimal_lexico_order(self, dictionary=None):
        closed_vertices = []
        to_add_to_dict = []
        for v in range(self.m):
            if self.is_closed(list_2_pow[v]):
                closed_vertices.append(v)
        minimal_facets_bin = self.facets_bin.copy()
        for v in closed_vertices:
            ridges_containing_v = [ridge for ridge in self.FP_bin[self.n - 2] if
                                   ridge | list_2_pow[v] == ridge]
            for ridge in ridges_containing_v:
                ridge_labels_to_modify = self.binary_to_face(ridge ^ list_2_pow[v])
                for ridge_labels_permu_iter in permutations(ridge_labels_to_modify):
                    F = []
                    for facet in self.facets_bin:
                        if facet | ridge == facet:
                            F += (self.binary_to_face(facet ^ ridge))
                    remaining_labels = self.binary_to_face(
                        (list_2_pow[self.m] - 1) ^ ridge ^ list_2_pow[F[0]] ^ list_2_pow[F[1]])
                    for F_perm_iter in permutations(F):
                        for remaining_labels_perm_iter in permutations(remaining_labels):
                            old_labels = [v]
                            old_labels += list(ridge_labels_permu_iter)
                            old_labels += list(F_perm_iter)
                            old_labels += list(remaining_labels_perm_iter)
                            relabeled_facets = relabel_facets(self, old_labels)
                            if dictionary != None:
                                if json.dumps(relabeled_facets) in dictionary:
                                    return relabeled_facets
                                else:
                                    to_add_to_dict.append(relabeled_facets)
                            if relabeled_facets < minimal_facets_bin:
                                minimal_facets_bin = relabeled_facets.copy()
        if dictionary != None:
            for relabeled_facets in to_add_to_dict:
                dictionary[json.dumps(relabeled_facets)] = False
            dictionary[json.dumps(minimal_facets_bin)] = True
        return minimal_facets_bin

    def is_a_seed(self):
        if self.MNF_set_bin == None:
            self.compute_MNF_set()
        list_all_edges = [face_to_binary(list(facet_iter), self.m) for facet_iter in
                          combinations(range(1, self.m + 1), 2)]
        for edge_bin in list_all_edges:
            is_pair = True
            for MNF_bin in self.MNF_set_bin:
                if MNF_bin & edge_bin in list_2_pow:  # 10 or 01 at the positions of the vertices of the edge
                    is_pair = False
            if is_pair and edge_bin not in self.MNF_set_bin:
                return False
        return True

    def give_non_seed_vertex(self):
        if self.MNF_set_bin == None:
            self.compute_MNF_set()
        list_all_edges = [face_to_binary(list(facet_iter), self.m) for facet_iter in
                          combinations(range(1, self.m + 1), 2)]
        for edge_bin in list_all_edges:
            is_pair = True
            for MNF_bin in self.MNF_set_bin:
                if MNF_bin & edge_bin in list_2_pow:  # 10 or 01 at the positions of the vertices of the edge
                    is_pair = False
            if is_pair and edge_bin not in self.MNF_set_bin:
                return binary_to_face(edge_bin, self.m)
        return []

    def relabel_canonical(self):
        if self.facets[0] != list(range(1, self.n + 1)):
            permutation_inv = (self.facets[0]).copy()
            for k in range(1, self.m + 1):
                if k not in permutation_inv:
                    permutation_inv.append(k)
            permutation = [0] * self.m
            for k in range(1, self.m + 1):
                permutation[permutation_inv[k - 1]-1] = k
            new_facets = []
            for facet in self.facets:
                new_facets.append([])
                for v in facet:
                    new_facets[-1].append(permutation[v - 1])
                new_facets[-1].sort()
            new_facets.sort()
            self.__init__(new_facets)


    def find_seed(self):
        return ()

def join(K,L):
    K.compute_MNF_set()
    L.compute_MNF_set()
    MNF_set_bin_join = K.MNF_set_bin.copy()
    for MNF in L.MNF_set_bin:
        MNF_set_bin_join.append(MNF<<K.m)
    return PureSimplicialComplex(None,[binary_to_face(MNF,K.m+L.m) for MNF in MNF_set_bin_join],K.n+L.n)

def wedge(K, v):
    if v > K.m or v < 1:
        raise Exception
    K.compute_MNF_set()
    K.MNF_bin_to_MNF()
    new_MNF_set = []
    for MNF in K.MNF_set:
        new_MNF_set.append([])
        for vertex in MNF:
            if vertex == v:
                new_MNF_set[-1].append(vertex)
                new_MNF_set[-1].append(K.m + 1)
            else:
                new_MNF_set[-1].append(vertex)
        new_MNF_set[-1].sort()
    new_MNF_set.sort()
    return PureSimplicialComplex(None, new_MNF_set, K.n + 1)


def multiple_wedge(K, J):
    new_K = PureSimplicialComplex(K.facets_bin)
    sum_wedge = 0
    if len(J) != K.m:
        raise Exception
    for v in range(K.m):
        for k in range(J[v]):
            new_K = wedge(new_K, v + 1)
            sum_wedge += 1
    new_K.relabel_canonical()
    return new_K


def give_IDCM_up_to_isom(K_J, J):
    all_IDCM = IDCM_Garrison_Scott(K_J)
    if len(all_IDCM) == 0:
        return []
    list_isom = [all_IDCM[0]]
    vertices_to_permute = []
    for r in range(K_J.m):
        if J[r] > 0:
            r_neighbors = (1 << K_J.m) - 1
            for MNF_bin in K_J.MNF_set_bin:
                if (1 << r) | MNF_bin == MNF_bin:
                    r_neighbors &= MNF_bin


def are_isom(K1:PureSimplicialComplex, K2:PureSimplicialComplex):
    if K1.n != K2.n or K1.m != K2.m or len(K1.facets_bin) != len(K2.facets_bin):
        return False
    # if (K1.facets_bin == K2.facets_bin).all:
    #     print('hello')
    #     return True
    K1.create_f_vector()
    K2.create_f_vector()
    if K1.f_vector != K2.f_vector:
        return False
    K1.compute_MNF_set()
    K2.compute_MNF_set()
    K1.MNF_bin_to_MNF()
    K2.MNF_bin_to_MNF()
    sizes_MNF_K1 = [len(MNF) for MNF in K1.MNF_set]
    sizes_MNF_K2 = [len(MNF) for MNF in K2.MNF_set]
    sizes_MNF_K1.sort()
    sizes_MNF_K2.sort()
    if sizes_MNF_K1 != sizes_MNF_K2:
        return False
    list_seq_K1 = []
    for vertex in range(1, K1.m + 1):
        list_seq_K1.append(sorted([len(MNF) for MNF in K1.MNF_set if vertex in MNF]))
    list_seq_K2 = []
    for vertex in range(1, K2.m + 1):
        list_seq_K2.append(sorted([len(MNF) for MNF in K2.MNF_set if vertex in MNF]))
    if sorted(list_seq_K1) != sorted(list_seq_K2):
        return False
    types_dict_K1 = dict()
    for index_vertex in range(K1.m):
        vertex = index_vertex + 1
        if json.dumps(list_seq_K1[index_vertex]) not in types_dict_K1:
            types_dict_K1[json.dumps(list_seq_K1[index_vertex])] = [vertex]
        else:
            types_dict_K1[json.dumps(list_seq_K1[index_vertex])].append(vertex)
    types_dict_K2 = dict()
    for index_vertex in range(K2.m):
        vertex = index_vertex + 1
        if json.dumps(list_seq_K2[index_vertex]) not in types_dict_K2:
            types_dict_K2[json.dumps(list_seq_K2[index_vertex])] = [vertex]
        else:
            types_dict_K2[json.dumps(list_seq_K2[index_vertex])].append(vertex)
    list_bij_K1 = []
    list_bij_K2 = []
    for item in types_dict_K1.items():
        seq, list_vertices = item
        list_bij_K1.append(list_vertices)
        list_bij_K2.append(types_dict_K2[seq])
    permutations_relabelling = []
    for k in range(len(list_bij_K2)):
        permutations_relabelling.append([list(permutation_iter) for permutation_iter in permutations(list_bij_K2[k])])
    K1_labels = []
    for labels in list_bij_K1:
        K1_labels += labels
    i = 0
    j = 0
    current_relabelling = []
    list_positions = [0] * len(permutations_relabelling)
    while i >= 0:
        if i == len(permutations_relabelling):
            K2_labels = []
            for labels in current_relabelling:
                K2_labels += labels
            help_bij = [(K1_labels[k], K2_labels[k]) for k in range(len(K1_labels))]
            help_bij.sort()
            old_labels = [data[1] - 1 for data in help_bij]
            if (K1.MNF_set_bin == relabel_MNF(K2, old_labels)).all():
                return True
            if i != 0:
                current_relabelling.pop()
            i -= 1
            list_positions[i] += 1
            continue
        j = list_positions[i]
        if j == len(permutations_relabelling[i]):
            list_positions[i] = 0
            if i != 0:
                current_relabelling.pop()
            i -= 1
            list_positions[i] += 1
            continue
        else:
            current_relabelling.append(permutations_relabelling[i][j])
            list_positions[i] = j
            i += 1
    return False


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


list_n_max = [2, 3, 5, 11]
seed_DB = []


for pic in range(1, 5):
    seed_DB.append([])
    if pic == 1:
        start = 1
    else:
        start = 2
    for n in range(start, list_n_max[pic - 1]):
        seed_DB[-1].append([])
        results_path = 'final_results/CSPLS_%d_%d' % (n, n + pic)
        list_m_n_seeds = [json.loads(facets_bytes) for facets_bytes in read_file(results_path)]
        for facets_seed in list_m_n_seeds:
            seed_DB[-1][-1].append(PureSimplicialComplex(facets_seed))


def is_PLS_new(K):
    list_of_links = []
    for v in list_2_pow[:K.m]:
        link_v_K = Link_of(K, v)
        list_non_seed_vertices = link_v_K.give_non_seed_vertex()
        while list_non_seed_vertices:
            link_v_K = Link_of(link_v_K, list_2_pow[list_non_seed_vertices[0] - 1])
            list_non_seed_vertices = link_v_K.give_non_seed_vertex()
        list_of_links.append(PureSimplicialComplex(link_v_K.facets_bin))
    is_PLS = False
    for link_K in list_of_links:
        is_PLS = False
        for L in seed_DB[link_K.Pic - 1][link_K.n - 2]:
            if are_isom(link_K, L):
                is_PLS = True
                break
        if not is_PLS:
            break
    del list_of_links
    # if is_PLS:
    #     K = list_of_links[0]
    #     K.compute_MNF_set()
    #     K.MNF_bin_to_MNF()
    #     print(K.MNF_set)
    return is_PLS


def relabel_facets(K, old_labels):
    if len(old_labels) != K.m:
        raise Exception
    new_facets = [0] * len(K.facets_bin)
    for k in range(len(new_facets)):
        for l in range(K.m):
            if K.facets_bin[k] | list_2_pow[old_labels[l]] == K.facets_bin[k]:
                # the case when 01 or 10
                new_facets[k] += list_2_pow[l]
    new_facets.sort()
    return new_facets


def relabel_MNF(K, old_labels):
    if len(old_labels) != K.m:
        raise Exception
    new_MNF = np.zeros(len(K.MNF_set_bin))
    for k in range(len(new_MNF)):
        for l in range(K.m):
            if K.MNF_set_bin[k] | list_2_pow[old_labels[l]] == K.MNF_set_bin[k]:
                # the case when 01 or 10
                new_MNF[k] += list_2_pow[l]
    new_MNF.sort()
    return new_MNF


def suspension(K):
    new_facets = []
    for facet in K.facets_bin:
        new_facets.append((facet << 1) | 1)
    return PureSimplicialComplex(new_facets)


def Link_of(K, F):
    if not K.FP_bin:
        K.create_FP()
    k = 0
    l = 0
    while F % list_2_pow[l] != F:
        if list_2_pow[l] | F == F:
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
    complementary_faces = [face for face in K.FP_bin[K.n - k - 1]]
    for complementary_face in complementary_faces:
        if (complementary_face | F) in K.facets_bin:
            facets_of_Link.append(complementary_face)
    unused_labels = []
    for i in range(K.m):
        unused = True
        for facet in facets_of_Link:
            if list_2_pow[i] | facet == facet:
                unused = False
                break
        if unused:
            unused_labels.append(i)
    for i in range(len(unused_labels) - 1, -1, -1):
        for j in range(len(facets_of_Link)):
            unshifted_bits = facets_of_Link[j] % list_2_pow[unused_labels[i]]
            facets_of_Link[j] = ((facets_of_Link[j] ^ unshifted_bits) >> 1) ^ unshifted_bits
    if len(unused_labels) == 0:
        return PureSimplicialComplex(facets_of_Link)
    m_of_link = K.m
    label_unused = True
    while m_of_link > 0 and label_unused:
        for facet in facets_of_Link:
            if facet ^ list_2_pow[m_of_link - 1] == facet:
                label_unused = False
                break
        m_of_link -= 1
    return PureSimplicialComplex(facets_of_Link)


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


def Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m):
    if len(pile) >= 1:
        K = pile.pop()
        if not K.is_closed():
            if len(K.facets_bin) < G_vector[aimed_m - 5]:
                for face in candidate_facets_ref:
                    if K.is_candidate(face, aimed_m):
                        new_K = PureSimplicialComplex(K.facets_bin + [face])
                        if new_K.is_promising():
                            pile.append(new_K)
                            Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)
                        else:
                            Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)
        else:
            if K.is_Z2_homology_sphere() and K.is_minimal_lexico_order():
                results.append(K)
            Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)


def Garrison_Scott(K):
    while K.facets[0] != list(range(1, K.n + 1)):
        vertex_to_relabel = 0
        image_of_vertex = 0
        for i in range(1, K.n + 1):
            if K.facets[0][i - 1] != i:
                vertex_to_relabel = K.facets[0][i - 1]
                image_of_vertex = i
                break
        new_facets = []
        for MF in K.facets:
            new_facets.append([])
            for vertex in MF:
                if vertex == vertex_to_relabel:
                    new_facets[-1].append(image_of_vertex)
                elif vertex == image_of_vertex:
                    new_facets[-1].append(vertex_to_relabel)
                else:
                    new_facets[-1].append(vertex)
            new_facets[-1].sort()
        new_facets.sort()
        K = PureSimplicialComplex(new_facets)
        print(K.facets[0])
    K.create_FP()
    list_char_funct = []
    '''we create all the non zero elements of Z_2^n'''
    list_2_pow = [1]
    for k in range(1, K.m + 1):
        list_2_pow.append(list_2_pow[-1] * 2)
    list_S = []
    for k in range(K.n, K.m):
        list_S.append(list(range(1, list_2_pow[K.n])))
    K_full_FP = []
    for k in range(K.n):
        for face in K.FP_bin[k]:
            K_full_FP.append(face)
    current_lambda = [0] * K.m
    i = K.n
    for k in range(K.n):
        current_lambda[k] = list_2_pow[k]
    while i >= K.n:
        list_S[K.n - i] = list(range(1, list_2_pow[K.n]))
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
                    if list_S[K.n - i] == []:
                        break
        while i >= K.n:
            if list_S[K.n - i] == []:
                i -= 1
                continue
            current_lambda[i] = list_S[K.n - i].pop()
            if i + 1 == K.m:
                list_char_funct.append(current_lambda.copy())
                continue
            else:
                i += 1
                break
    return list_char_funct


def IDCM_Garrison_Scott(K):
    K.create_FP()
    m = K.m
    n = K.n
    p = K.Pic
    full_face = list_2_pow[m] - 1
    max_faces = [face ^ full_face for face in K.facets_bin]
    # We create a simplicial complex having its maximal faces represented by the cofacets of K
    L = PureSimplicialComplex(max_faces)
    L.create_FP()
    CF_full_set = []
    for k in range(L.n):
        for face_bin in L.FP_bin[k]:
            CF_full_set.append(face_bin)
    # The list CF_full_set represents all the sub-cofacets of K
    current_IDCM = [0] * m
    for k in range(n, m):
        current_IDCM[k] = list_2_pow[k - n]
    S = [k for k in range(1, list_2_pow[p])]
    list_S = [S.copy() for k in range(m)]
    list_results = []
    list_condition_i = [0] * m
    sum = 0
    for k in range(m - 1, -1, -1):
        sum += list_2_pow[k]
        list_condition_i[k] = sum
    i = n - 1
    while i < n:
        list_S[i] = S.copy()
        for c in current_IDCM[i + 1:]:
            list_S[i].remove(c)
        for CF in CF_full_set:
            if CF | list_condition_i[i] == list_condition_i[i] and CF | list_2_pow[i] == CF:
                list_indexes_CF = [j for j in range(m) if list_2_pow[j] | CF == CF and j != i]
                linear_comb = 0
                for j in list_indexes_CF:
                    linear_comb ^= current_IDCM[j]
                if linear_comb in list_S[i]:
                    list_S[i].remove(linear_comb)
                    if list_S[i] == []:
                        break
        while i < n:
            if list_S[i] == []:
                i += 1
                continue
            current_IDCM[i] = list_S[i].pop()
            if i == 0:
                list_results.append(current_IDCM.copy())
                continue
            else:
                i -= 1
                break
    return list_results


def int_to_bin_array(x, k):
    res = np.zeros(k)
    for i in range(k):
        if list_2_pow[i] | x == x:
            res[i] = 1
    return res


def char_funct_array_to_bin(char_funct):
    m, Pic = char_funct.shape
    np_list_2_pow = np.array(list_2_pow[:Pic])
    res = [np.sum(np_list_2_pow[v == 1]) for v in char_funct]
    return res


def enumerate_char_funct_orbits(n, m):
    list_char_funct_bin = []
    list_char_funct = []
    list_of_rows = []
    for k in range(3, 2 ** (m - n)):
        if k not in list_2_pow:
            list_of_rows.append(k)
    for combi_iter in combinations(list_of_rows, n):
        char_funct = list(combi_iter) + list_2_pow[:m - n]
        list_char_funct_bin.append(char_funct)
        current_char_funct = np.zeros((m, m - n))
        for i in range(m):
            current_char_funct[i] = int_to_bin_array(char_funct[i], m - n)
        list_char_funct.append(current_char_funct.copy())
    # GL4 = enumerate_GL4()
    SLn = enumerate_SL(m - n)
    eq_classes_repres = [list_char_funct_bin[0]]
    eq_classes_ref = []
    mini = sorted(list_char_funct_bin[0])
    for A in SLn:
        new_char_funct = list_char_funct[0].dot(A)
        new_char_funct %= 2
        min_new_char_funct = sorted(char_funct_array_to_bin(new_char_funct))
        if min_new_char_funct < mini:
            mini = min_new_char_funct.copy()
    eq_classes_ref.append(mini)
    for i in tqdm(range(1, len(list_char_funct))):
        is_a_new_repres = True
        mini = sorted(list_char_funct_bin[i])
        for A in SLn:
            new_char_funct = list_char_funct[i].dot(A)
            new_char_funct %= 2
            min_new_char_funct = sorted(char_funct_array_to_bin(new_char_funct))
            # keeps_colours = True
            # for element in list_2_pow[:m-n-3]:
            #     if element not in min_new_char_funct:
            #         keeps_colours = False
            #         break
            # if not keeps_colours:
            #     break
            if min_new_char_funct < mini:
                mini = min_new_char_funct.copy()
            if min_new_char_funct in eq_classes_ref:
                is_a_new_repres = False
                break
        if is_a_new_repres:
            eq_classes_repres.append(list_char_funct_bin[i])
            eq_classes_ref.append(mini)
    return (eq_classes_repres)

def find_facets_compatible_with_lambda(char_function, m, n):
    Pic = m - n
    cofacets = []
    for cofacet_iter in combinations(range(1, m + 1), Pic):
        sub_array = []
        for index in cofacet_iter:
            sub_array.append(char_function[index - 1])
        if Z2la.Z2Array(Pic, sub_array.copy()).is_invertible():
            cofacets.append(list(cofacet_iter))
    facets = []
    for cofacet in cofacets:
        facets.append((list_2_pow[m] - 1) ^ face_to_binary(cofacet, m))
    facets.sort()
    return facets

def enumerate_non_isom_bin_matroids(n,m): #Computes the binary matroids up to isomorphisms
    list_lambdas = enumerate_char_funct_orbits(n,m)
    print(len(list_lambdas))
    list_M = [PureSimplicialComplex(find_facets_compatible_with_lambda(char_map,m,n)) for char_map in list_lambdas]
    list_nonisom_matroids = []
    for M1 in list_M:
        is_isom = False
        for M2 in list_nonisom_matroids:
            if are_isom(M1,M2):
                is_isom = True
                break
        if not is_isom:
            list_nonisom_matroids.append(M1)
            # M1.compute_MNF_set()
            # M1.MNF_bin_to_MNF()
            # print(M1.MNF_set)
    list_facets_nonisom_matroids = [M.facets_bin for M in list_nonisom_matroids]
    return list_facets_nonisom_matroids

def text(result):
    name = '/GL4'
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


def enumerate_GL4():
    GL4 = []
    for i0 in range(1, 16):
        for i1 in range(1, 16):
            for i2 in range(1, 16):
                for i3 in range(1, 16):
                    A_bin = Z2la.Z2Array(4, [i0, i1, i2, i3])
                    if A_bin.is_invertible():
                        A = np.zeros((4, 4))
                        A[0] = int_to_bin_array(i0, 4)
                        A[1] = int_to_bin_array(i1, 4)
                        A[2] = int_to_bin_array(i2, 4)
                        A[3] = int_to_bin_array(i3, 4)
                        GL4.append(A.copy())
    return (GL4)


def enumerate_SL(n):
    A = np.eye(n)
    SLn = []
    list_of_vectors = []
    for vect in A:
        list_of_vectors.append(vect)
    for vect in permutations(
            list_of_vectors):
        current_element = np.zeros((n, n))
        for k in range(n):
            current_element[k] = vect[k]
        SLn.append(current_element.copy())
    return SLn


def char_funct_bin_to_numpy(char_funct, n):
    char_funct_array = np.zeros((n, len(char_funct)))
    for k in range(len(char_funct)):
        char_funct_array[:, k] = int_to_bin_array(char_funct[k], n)
    return(char_funct_array)


def give_next_vect(vect, base):
    index = 0
    vect[index] = (vect[index] + 1) % base[index]
    while index < vect.size - 1 and vect[index] == 0:
        index += 1
        vect[index] = (vect[index] + 1) % base[index]


def find_Z4_homology(K, IDCM):
    n = K.n
    m = K.m
    H = np.zeros(K.n + 1)
    list_gene_omega = np.array([(IDCM[d] << n) + list_2_pow[d] for d in range(n)], dtype=int)
    vect = np.zeros(n)
    list_omega = np.zeros(2 ** n, dtype=int)
    for k in range(2 ** n):
        list_omega[k] = np.bitwise_xor.reduce(list_gene_omega[np.flatnonzero(vect)])
        give_next_vect(vect, 2 * np.ones(n))
    list_H_omega = []
    for omega in list_omega:
        if omega == 0:
            list_H_omega.append([])
        else:
            K_omega = PureSimplicialComplex([omega & facet_bin for facet_bin in K.facets_bin if omega & facet_bin != 0])
            K_omega.Z2_Betti_numbers()
            list_H_omega.append(K_omega.H)
    for H_omega in list_H_omega:
        # print(H_omega)
        if not H_omega:
            H[0] += 1
        else:
            for i in range(len(H_omega)):
                if i == 0:
                    H[i + 1] += (H_omega[i] - 1)
                else:
                    H[i + 1] += H_omega[i]
    return H

def DCM_bin_to_IDCM_bin(DCM_bin,n):
    m = len(DCM_bin)
    M = np.eye(m,dtype=int)
    for i in range(n+1,m):
        for j in range(m):
            if list_2_pow[j]|DCM_bin[i]==DCM_bin[i]:
                M[j,i]=1
    return(np.dot(M[:,n:],np.array((list_2_pow[:m-n]))))


def enumerate_chordless_cycles(A,i,j,list_cycles):
    def enumerate_chordless_cycles_rec(A, i0, j0, i, j, current_cycle, list_cycles, horizontal):
        if horizontal:
            for next_j in np.where(np.logical_and((current_cycle==0).all(axis=0),(A-current_cycle)[i,:]==1))[0]:
                next_cycle = np.copy(current_cycle)
                next_cycle[i, next_j] = 1
                enumerate_chordless_cycles_rec(A, i0, j0, i, next_j, next_cycle, list_cycles, False)
        else:
            if j==j0 and (A-current_cycle)[i0,j]==1:
                current_cycle[i0,j0]=1
                filter_i = (current_cycle==1).any(axis=1)
                filter_j = (current_cycle==1).any(axis=0)
                if (current_cycle[filter_i,:][:,filter_j]==A[filter_i,:][:,filter_j]).all():
                    list_cycles.append(current_cycle)
            else:
                for next_i in np.where(np.logical_and((current_cycle==0).all(axis=1),(A-current_cycle)[:,j]==1))[0]:
                    next_cycle = np.copy(current_cycle)
                    next_cycle[next_i,j] = 1
                    enumerate_chordless_cycles_rec(A, i0, j0, next_i, j, next_cycle, list_cycles, True)
    n,m = A.shape
    enumerate_chordless_cycles_rec(A,i,j,i,j,np.zeros((n,m)),list_cycles,True)

def lifting_algo(original_facets,char_map):
    # INITIALIZATION #
    n,m = char_map.shape
    relabeled_facets = np.zeros((len(original_facets),m),dtype=bool)
    facets = [original_facets[0]]
    appeared = np.ones(m+1,dtype = bool)
    appeared[:n+1] = False
    unused_facets = []
    for k in range(1,len(original_facets)):
        facet = original_facets[k]
        if np.count_nonzero(appeared[facet])<=1:
            facets.append(facet)
            appeared[facet] = False
        else:
            unused_facets.append(facet)
    for facet in unused_facets:
        facets.append(facet)
    for k in range(len(facets)):
        facet = facets[k]
        for v in facet:
            relabeled_facets[k][v-1] = True
    appeared = [False]*m
    ordered_vertices = []
    for facet in relabeled_facets:
        for v in np.where(facet)[0]:
            if not appeared[v]:
                appeared[v] = True
                ordered_vertices.append(v)
    ordered_vertices = np.array(ordered_vertices)
    # END OF INITIALIZATION #
    tree = np.zeros((n,m))
    tree[:,:n+1] = char_map[:,:n+1]
    for k in range(n+1,m):
        tree[(tree[:,n:k]==0).all(axis=1),k] = char_map[(tree[:,n:k]==0).all(axis=1),k]
        rows_having_an_edge_previously = (tree[:, n:k] == 1).any(axis=1)
        rows_having_an_edge_now = (char_map[:,k]==1)
        if np.logical_and(rows_having_an_edge_previously,rows_having_an_edge_now).any():
            index = np.where(np.logical_and(rows_having_an_edge_previously,rows_having_an_edge_now))[0][0]
            tree[index,k]=1
    undecided_edges = char_map - tree
    relabelling = char_map - undecided_edges
    stack_cycles_conflict = []
    for facet in relabeled_facets:
        facet_subgraph = np.zeros((n,m))
        facet_subgraph[:,facet] = char_map[:,facet]
        undecided_edges_subgraph = np.zeros((n,m))
        undecided_edges_subgraph[:, facet] = undecided_edges[:,facet]
        positions = np.where(undecided_edges_subgraph)
        for k in range(len(positions[0])):
            i = positions[0][k]
            j = positions[1][k]
            list_cycles = []
            enumerate_chordless_cycles(facet_subgraph,i,j,list_cycles)
            for cycle in list_cycles:
                if np.count_nonzero(np.multiply(cycle,undecided_edges_subgraph))==1:
                    if relabelling[i,j]==0:
                        print("completing cycle:")
                        print(cycle)
                        if np.sum(np.multiply(cycle,relabelling))%4==1:
                            relabelling[i,j] = -1
                        elif np.sum(np.multiply(cycle,relabelling))%4==3:
                            relabelling[i,j] = 1
                        else:
                            print("bizarre")
                    elif np.sum(np.multiply(cycle,relabelling))%4==2:
                        print("Conflict middle!")
                else:
                    stack_cycles_conflict.append(cycle)
    if (np.abs(relabelling)==char_map).all():
        conflict = False
        for cycle in stack_cycles_conflict:
            if np.count_nonzero(np.multiply(cycle, relabelling))%4==2:
                conflict = True
                print(cycle)
                break
        if conflict:
            print("Conflict end!")
            print(undecided_edges)
            print(char_map)
            print(facets)

    return(relabelling)

def integer_det(M):
    M = [row[:] for row in M] # make a copy to keep original M unmodified
    N, sign, prev = len(M), 1, 1
    for i in range(N-1):
        if M[i][i] == 0: # swap with another row having nonzero i's elem
            swapto = next( (j for j in range(i+1,N) if M[j][i] != 0), None )
            if swapto is None:
                return 0 # all M[*][i] are zero => zero determinant
            M[i], M[swapto], sign = M[swapto], M[i], -sign
        for j in range(i+1,N):
            for k in range(i+1,N):
                assert ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) % prev == 0
                M[j][k] = ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) // prev
        prev = M[i][i]
    return sign * M[-1][-1]
def integer_det_mod_2(M):
    M = [row[:] for row in M] # make a copy to keep original M unmodified
    N, sign, prev = len(M), 1, 1
    for i in range(N-1):
        if M[i][i] == 0: # swap with another row having nonzero i's elem
            swapto = next( (j for j in range(i+1,N) if M[j][i] != 0), None )
            if swapto is None:
                return 0 # all M[*][i] are zero => zero determinant
            M[i], M[swapto], sign = M[swapto], M[i], -sign
        for j in range(i+1,N):
            for k in range(i+1,N):
                M[j][k] = (( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) * prev) % 2
        prev = M[i][i]
    return M[-1][-1]


def lifting_binary_matroid(M):
    n,m = M.shape
    bases = []
    for iter in combinations(range(m),n):
        if integer_det_mod_2(M[:, iter]) == 1:
            bases.append(list(iter))
    tree = np.zeros((n,m))
    tree[:,:n] = M[:,:n]
    for k in range(n,m):
        tree[(tree[:,n:k]==0).all(axis=1),k] = M[(tree[:,n:k]==0).all(axis=1),k]
        rows_having_an_edge_previously = (tree[:, n:k] == 1).any(axis=1)
        rows_having_an_edge_now = (M[:,k]==1)
        if np.logical_and(rows_having_an_edge_previously,rows_having_an_edge_now).any():
            index = np.where(np.logical_and(rows_having_an_edge_previously,rows_having_an_edge_now))[0][0]
            tree[index,k]=1
    undecided_edges = M - tree
    positions = np.where(undecided_edges)
    list_positions=[(positions[0][k],positions[1][k]) for k in range(len(positions[0]))]
    N = len(list_positions)
    for k in range(0,N+1):
        for iter in combinations(list_positions,k):
            N = M.copy()
            for pair in iter:
                N[pair[0],pair[1]] = -1
            is_a_lift = True
            list_problematic_bases=[]
            list_dets=[]
            for base in bases:
                d=integer_det(N[:,base])
                if abs(d)!=1:
                    is_a_lift= False
                    list_problematic_bases.append(base)
                    list_dets.append(d)
                    # break
            if len(list_problematic_bases)<=5:
                print(np.array(range(m)))
                print(N)
                print(list_problematic_bases)
                print(list_dets)
            if is_a_lift:
                return True
    return False



