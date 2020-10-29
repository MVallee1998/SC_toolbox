import sqlite3
import numpy.polynomial as npp
import binary_search_tree as bst
import json
import datetime
from itertools import combinations, permutations
import timeit


def face_to_binary(face, m):
    binary_face = m * [0]
    for k in face:
        binary_face[m - k] = 1
    return int(''.join(map(str, binary_face)), 2)


class PureSimplicialComplex:
    def __init__(self, n, m, facets=None, MNF_set=None):
        self.n = n
        self.m = m
        self.MNF_set = MNF_set
        self.facets = facets
        if len(facets) > 0:
            self.facets_bin = [face_to_binary(facet, self.m) for facet in self.facets]
            self.faces_set_to_MNF_set()
        elif len(MNF_set) > 0:
            self.MNF_set_bin = [face_to_binary(MNF, self.m) for MNF in self.MNF_set]
            self.compute_facets_from_MNF_set()
        self.FP_bin = None
        self.Pic = m - n
        self.f_vector = None
        self.h_vector = None
        self.g_vector = None
        self.list_2_pow = [1]
        for k in range(m + 1):
            self.list_2_pow.append(self.list_2_pow[-1] * 2)

    def create_FP(self):
        self.FP_bin = [[]] * self.n
        faces_set = [bst.Node() for k in range(self.n)]
        facet_data = [(facet, []) for facet in self.facets]
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
                        if not faces_set[k].findval((subface, [l])):
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
        # we enumerate all the faces until the time we will have the right dimension faces
        while k > self.n:
            subface_set = bst.Node()
            for face in temp_faces_list:  # treat every face of dimension k
                for element in self.list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face == face:
                        subface = element ^ face
                        if not subface_set.findval(subface):
                            if k > self.n or (k == self.n and (not dichotomie(subface, self.FP_bin[k - 1]))):
                                subface_set.insert(subface)
            temp_faces_list = []
            subface_set.TreeToList(temp_faces_list)
            k = k - 1
        # we then create the set with all subsets of [m] of size <= n not in  the given faces set
        all_non_faces = [bst.Node() for k in range(self.n + 1)]
        all_non_faces[-1] = subface_set
        for k in range(self.n - 1, -1, -1):  # dimension
            faces = []
            all_non_faces[k + 1].TreeToList(faces)
            for face in faces:  # treat every face of dimension k
                for element in self.list_2_pow[:self.m]:  # construct the (k-1)-subfaces
                    if element | face == face:
                        subface = element ^ face
                        if (not (all_non_faces[k].findval(subface))) and (not dichotomie(subface, self.FP_bin[k])):
                            all_non_faces[k].insert(subface)
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
                    MNF_set.insert(face_to_test)
        self.MNF_set_bin = []
        MNF_set.TreeToList(self.MNF_set_bin)

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

    def is_closed(self):
        if not self.FP_bin:
            self.create_FP()
        if self.n < 2: return True
        is_closed = True
        for ridge_data in self.FP_bin[self.n - 2]:
            if len(ridge_data[1]) != 2:
                is_closed = False
                break
        return is_closed

    def is_minimal_lexico_order(self):
        closed_vertices = []
        for v in range(self.m):
            link_of_v = Link_of(self, self.list_2_pow[v])
            if link_of_v.is_closed():
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
    facets_of_K = [facet_data[0] for facet_data in K.FP_bin[K.n - 1]]
    k = len(F) - 1
    is_a_face = False
    for facet in facets_of_K:
        if F | facet == facet:
            is_a_face = True
            break
    if not is_a_face:
        raise Exception
    facets_of_Link = []
    complementary_faces = [facet_data[0] for facet_data in K.FP_bin[K.n - k - 1]]
    for complementary_face in complementary_faces:
        if dichotomie(complementary_face ^ F, facets_of_K):
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
        return PureSimplicialComplex(K.n - k - 1, K.m, facets_of_Link)
    m_of_link = K.m
    label_unused = True
    while m_of_link > 0 & label_unused:
        for facet in facets_of_Link:
            if facet ^ K.list_2_pow[m_of_link - 1] == facet:
                label_unused = False
                break
    return PureSimplicialComplex(K.n - k - 1, m_of_link, facets_of_Link)


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


def Garrison_Scott(K=PureSimplicialComplex()):
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
