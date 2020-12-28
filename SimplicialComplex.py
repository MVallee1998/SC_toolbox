import numpy.polynomial as npp
import binary_search_tree as bst
from itertools import combinations, permutations
import Betti_numbers as bnbr
import numpy as np
import sys
import Z2_linear_algebra as Z2la
import json

# sys.setrecursionlimit(1000)

list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1] * 2)


def face_to_binary(face, m):
    binary_face = m * [0]
    for k in face:
        binary_face[m - k] = 1
    return int(''.join(map(str, binary_face)), 2)


def binary_to_face(x, m):
    result = []
    for k in range(m):
        if list_2_pow[k] | x == x:
            result.append(k + 1)
    return result


class PureSimplicialComplex:
    def __init__(self, facets=None, MNF_set=None, n=0):
        self.m = 0
        self.n = n
        self.facets = facets.copy()
        self.MNF_set = MNF_set
        self.facets_bin = []
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
                self.facets_bin = [face_to_binary(facet, self.m) for facet in self.facets]
            else:
                self.n = 0
                self.m = 0
                self.facets_bin = facets.copy()
                self.facets_bin.sort()
                last_facet = self.facets_bin[-1]
                two_pow = 1
                while last_facet % two_pow != last_facet:
                    if two_pow | last_facet == last_facet:
                        self.n += 1
                    two_pow *= 2
                for element in list_2_pow:
                    label_used = False
                    for facet in self.facets_bin:
                        if facet | element == facet:
                            label_used = True
                            break
                    if label_used:
                        self.m += 1
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
        self.FP_bin = None
        self.Pic = self.m - self.n
        self.f_vector = None
        self.h_vector = None
        self.g_vector = None

        self.H = None
        self.unclosed_ridges = None
        self.closed_faces = None

    def create_FP(self):
        if self.facets_bin and not self.FP_bin:
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
                    for element in list_2_pow[:self.m]:  # construct the (k-1)-subfaces
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
        if not self.FP_bin:
            self.create_FP()
        temp_faces_list = [list_2_pow[self.m] - 1]
        k = self.m
        subface_set = bst.Node()
        subface_set.insertList(temp_faces_list)
        # we enumerate all the faces until the time we will have the right dimension faces
        while k > self.n:
            subface_set = bst.Node()
            for face in temp_faces_list:  # treat every face of dimension k
                for element in list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face == face:
                        subface = element ^ face
                        if not subface_set.findval(subface):
                            if k > self.n or (k == self.n and (not dichotomie(self.faces_set[k - 1], subface))):
                                subface_set.insert((subface, []))
            temp_faces_list_data = []
            subface_set.TreeToList(temp_faces_list_data)
            temp_faces_list = [data[0] for data in temp_faces_list_data]

            k = k - 1
        # we then create the set with all subsets of [m] of size <= n not in  the given faces set
        all_non_faces = [bst.Node() for k in range(self.n + 1)]
        all_non_faces[-1] = subface_set
        for k in range(self.n - 1, -1, -1):  # dimension
            faces_data = []
            all_non_faces[k + 1].TreeToList(faces_data)
            faces = [data[0] for data in faces_data]
            for l in range(len(faces)):  # treat every face of dimension k
                face = faces[l]
                for element in list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face == face:
                        subface = element ^ face
                        if (not (all_non_faces[k].findval(subface))) and (
                                not dichotomie(self.FP_bin[k][:][0], subface)):  # needs to be fixed
                            all_non_faces[k].insert((subface, [l + 1]))
        MNF_set = bst.Node()
        # Here we will try every face of dimension k
        for k in range(1, self.n + 1):
            faces_to_test_data = []
            all_non_faces[k].TreeToList(faces_to_test_data)
            faces_to_test = [data[0] for data in faces_to_test_data]
            for face_to_test in faces_to_test:
                is_MNF = True
                for element in list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face_to_test == face_to_test:
                        subface = element ^ face_to_test
                        if not (dichotomie(self.FP_bin[k - 1][:][0], subface)):
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
            candidate_facets.append(sum([list_2_pow[k] for k in facet_iter]))
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
            if face_bin | list_2_pow[k] == face_bin:
                face.append(k)
        return face

    def is_closed(self, F=None):
        if F:
            Link_of_F = Link_of(self, F)
            ridges_set = bst.Node()
            if Link_of_F.n < 2: return True
            for k in range(len(Link_of_F.facets_bin)):
                facet = Link_of_F.facets_bin[k]
                for element in list_2_pow:
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
                    for element in list_2_pow:
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
            for k in range(self.n - 1):
                for face_data in self.FP_bin[k]:
                    if self.is_closed(face_data[0]):
                        self.closed_faces.append(face_data[0])

    def list_unclosed_ridges(self):
        if not self.unclosed_ridges:
            ridges_set = bst.Node()
            if self.n < 2: return []
            for k in range(len(self.facets_bin)):
                facet = self.facets_bin[k]
                for element in list_2_pow:
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
        for element in list_2_pow:
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
            if not self.facets_bin:
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
            if self.is_closed(list_2_pow[v]):
                closed_vertices.append(v)
        minimal_facets_bin = self.facets_bin.copy()
        for v in closed_vertices:
            ridges_containing_v = [ridge_data[0] for ridge_data in self.FP_bin[self.n - 2] if
                                   ridge_data[0] | list_2_pow[v] == ridge_data[0]]
            for ridge in ridges_containing_v:
                ridge_labels_to_modify = self.binary_to_face(ridge ^ list_2_pow[v])
                for ridge_labels_permu_iter in permutations(ridge_labels_to_modify):
                    F = []
                    for facet in self.facets_bin:
                        if facet | ridge == facet:
                            F += (self.binary_to_face(facet ^ ridge))
                    remaining_labels = self.binary_to_face((list_2_pow[self.m] - 1) ^ ridge ^ list_2_pow[F[0]] ^ list_2_pow[F[1]])
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

    def find_minimal_lexico_order(self,dictionary):
        closed_vertices = []
        for v in range(self.m):
            if self.is_closed(list_2_pow[v]):
                closed_vertices.append(v)
        minimal_facets_bin = self.facets_bin.copy()
        for v in closed_vertices:
            ridges_containing_v = [ridge_data[0] for ridge_data in self.FP_bin[self.n - 2] if
                                   ridge_data[0] | list_2_pow[v] == ridge_data[0]]
            for ridge in ridges_containing_v:
                ridge_labels_to_modify = self.binary_to_face(ridge ^ list_2_pow[v])
                for ridge_labels_permu_iter in permutations(ridge_labels_to_modify):
                    F = []
                    for facet in self.facets_bin:
                        if facet | ridge == facet:
                            F += (self.binary_to_face(facet ^ ridge))
                    remaining_labels = self.binary_to_face((list_2_pow[self.m] - 1) ^ ridge ^ list_2_pow[F[0]] ^ list_2_pow[F[1]])
                    for F_perm_iter in permutations(F):
                        for remaining_labels_perm_iter in permutations(remaining_labels):
                            old_labels = [v]
                            old_labels += list(ridge_labels_permu_iter)
                            old_labels += list(F_perm_iter)
                            old_labels += list(remaining_labels_perm_iter)
                            relabeled_facets = relabel_facets(self, old_labels)
                            if json.dumps(relabeled_facets) in dictionary:
                                dictionary[json.dumps(relabeled_facets)] = True
                            if relabeled_facets < minimal_facets_bin:
                                minimal_facets_bin = relabeled_facets.copy()
        return minimal_facets_bin




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
    complementary_faces = [facet_data[0] for facet_data in K.FP_bin[K.n - k - 1]]
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
    if pile == []:
        return results
    K = pile.pop()
    if not K.is_closed():
        for face in candidate_facets_ref:
            new_K = PureSimplicialComplex(K.facets_bin + [face])
            pile.append(new_K)
        return Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)
    else:
        results.append(K)
        print(K.facets_bin)
        return Hyuntae_algo(pile, candidate_facets_ref, results, aimed_m)


def Garrison_Scott(K):
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
            K_full_FP.append(face[0])
    current_lambda = []
    i = K.n
    for k in range(K.n):
        current_lambda.append(list_2_pow[k])
    for k in range(K.n, K.m):
        current_lambda.append(0)
    while i >= K.n:
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


def int_to_bin_array(x, k):
    res = np.zeros(k)
    for i in range(k):
        if list_2_pow[i] | x == x:
            res[i] = 1
    return res


def char_funct_array_to_bin(char_funct):
    m, Pic = char_funct.shape
    res = []
    np_list_2_pow = np.array(list_2_pow[:Pic])
    for v in char_funct:
        res.append(np.sum(np_list_2_pow[v == 1]))
    return res


def enumerate_char_funct_orbits(n, m):
    list_char_funct_bin = []
    list_char_funct = []
    for combi_iter in combinations([3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15], n):
        char_funct = list(combi_iter) + list_2_pow[:m - n]
        list_char_funct_bin.append(char_funct)
        current_char_funct = np.zeros((m, m - n))
        for i in range(m):
            current_char_funct[i] = int_to_bin_array(char_funct[i], m - n)
        list_char_funct.append(current_char_funct.copy())
    # GL4 = enumerate_GL4()
    SL4 = enumerate_SL4()
    eq_classes_repres = [list_char_funct_bin[0]]
    eq_classes_ref = []
    mini = sorted(list_char_funct_bin[0])
    for A in SL4:
        new_char_funct = list_char_funct[0].dot(A)
        new_char_funct %= 2
        min_new_char_funct = sorted(char_funct_array_to_bin(new_char_funct))
        if min_new_char_funct < mini:
            mini = min_new_char_funct.copy()
    eq_classes_ref.append(mini)
    for i in range(1, len(list_char_funct)):
        is_a_new_repres = True
        mini = sorted(list_char_funct_bin[i])
        for A in SL4:
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

def enumerate_all_lambdas(n,m):
    list_char_funct = []
    for combi_iter in combinations([3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15], n):
        char_funct = list(combi_iter) + list_2_pow[:m - n]
        list_char_funct.append(char_funct)
    return list_char_funct

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


def enumerate_SL4():
    SL4 = []
    for vect in permutations(
            [np.array([1, 0, 0, 0]), np.array([0, 1, 0, 0]), np.array([0, 0, 1, 0]), np.array([0, 0, 0, 1])]):
        current_element = np.zeros((4, 4))
        for k in range(4):
            current_element[k] = vect[k]
        SL4.append(current_element.copy())
    return SL4

def display_char_funct(char_funct,n):
    char_funct_array = np.zeros((n,len(char_funct)))
    for k in range(len(char_funct)):
        char_funct_array[:,k] = int_to_bin_array(char_funct[k],n)
    print(char_funct_array)

# K = [31, 47, 55, 59, 93, 94, 109, 117, 174, 185, 186, 206, 217, 218, 233, 236, 310, 342, 369, 372, 422, 433, 434, 454, 465, 466, 481, 484]
# K_sp = PureSimplicialComplex(K)
#
# # print(K)
# print(K_sp.is_Z2_homology_sphere(),K_sp.is_promising(),K_sp.is_closed(),K_sp.is_minimal_lexico_order())
#
# print(Garrison_Scott(K_sp))
