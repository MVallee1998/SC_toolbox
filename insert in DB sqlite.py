import sqlite3
import numpy.polynomial as npp
import binary_search_tree as bst
import json
import datetime
from itertools import combinations, permutations


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


def create_faces_set(facets):
    n = len(facets[0])
    faces_set = [bst.Node() for k in range(n)]
    faces_set[-1].insertList(facets)
    for k in range(n - 2, -1, -1):  # dimension
        faces = []
        faces_set[k + 1].TreeToList(faces)
        for face in faces:  # treat every face of dimension k
            for i in range(len(face)):  # construct the (k-1)-subfaces
                subfacet = []
                for j in range(len(face)):  # for constructing the subfacet
                    if i != j:
                        subfacet.append(face[j])
                if not faces_set[k].findval(subfacet):
                    faces_set[k].insert(subfacet)
    return faces_set


def create_f_vector(faces_set):
    n = len(faces_set)
    f_vector = [1]
    for i in range(n):
        f_vector.append(faces_set[i].size())
    return f_vector


def create_h_vector(f_vector):
    f_vector_copy = list(f_vector)
    f_vector_copy.reverse()
    F_pol = npp.Polynomial(f_vector_copy)
    P = npp.Polynomial([-1, 1])
    H_pol = F_pol(P)
    h_vector = list(H_pol)
    h_vector = [int(h_i) for h_i in h_vector]
    return h_vector


def create_g_vector(h_vector):
    n = len(h_vector)
    g_vector = []
    for i in range(1, int(n / 2) + 1):
        g_vector.append(h_vector[i] - h_vector[i - 1])
    return g_vector


def faces_set_to_MNF_set(faces_set):
    facets = []
    faces_set[-1].TreeToList(facets)
    m = max([facet[-1] for facet in facets])
    n = len(facets[0])
    temp_faces_list = [list(range(1, m + 2))]
    k = m
    subfacets_set = bst.Node()
    # we enumerate all the faces until the time we will have the right dimension faces
    while k > n:
        subfacets_set = bst.Node()
        for face in temp_faces_list:  # treat every face of dimension k
            for i in range(len(face)):  # construct the (k-1)-subfaces
                subfacet = []
                for j in range(len(face)):  # for constructing the subfacet
                    if i != j:
                        subfacet.append(face[j])
                if not subfacets_set.findval(subfacet):
                    if k > n or (k == n and (not faces_set[k - 1].findval(subfacet))):
                        subfacets_set.insert(subfacet)
        temp_faces_list = []
        subfacets_set.TreeToList(temp_faces_list)
        k = k - 1
    # we then create the set with all subsets of [m] of size <= n not in  the given faces set
    all_non_faces = [bst.Node() for k in range(n + 1)]
    all_non_faces[-1] = subfacets_set
    for k in range(n - 1, -1, -1):  # dimension
        faces = []
        all_non_faces[k + 1].TreeToList(faces)
        for face in faces:  # treat every face of dimension k
            for i in range(len(face)):  # construct the (k-1)-subfaces
                subfacet = []
                for j in range(len(face)):  # for constructing the subfacet
                    if i != j:
                        subfacet.append(face[j])
                if (not (all_non_faces[k].findval(subfacet))) and (not faces_set[k].findval(subfacet)):
                    all_non_faces[k].insert(subfacet)
    MNF_set = bst.Node()
    # Here we will try every face of dimension
    for k in range(1, n + 1):
        faces_to_test = []
        all_non_faces[k].TreeToList(faces_to_test)
        for face_to_test in faces_to_test:
            is_MNF = True
            for i in range(len(face_to_test)):  # construct the (k-1)-subfaces
                subfacet = []
                for j in range(len(face_to_test)):  # for constructing the subfacet
                    if i != j:
                        subfacet.append(face_to_test[j])
                if not (faces_set[k - 1].findval(subfacet)):
                    is_MNF = False
                    break
            if is_MNF:
                MNF_set.insert(face_to_test)
    return MNF_set


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



def MNF_set_to_Maximal_faces_set(MNF_set, n, m):
    M = range(1, m + 1)
    candidate_facets_iter = combinations(M, n)
    facets = []
    for facet_iter in candidate_facets_iter:
        facet = list(facet_iter)
        is_a_facet = True
        for MNF in MNF_set:  # We chack if it the facet does not contain some minimal non-face
            if not len(MNF) > len(facet):
                is_contained = True
                for i in MNF:
                    if not dicho_search(i, facet, 0, len(facet) - 1):
                        is_contained = False
                        break
                if is_contained:
                    is_a_facet = False
        if is_a_facet:
            facets.append(facet)
    return facets




def find_minimal_facets_set(facets_set):
    facets_set.sort()
    n = len(facets_set[0])
    # if the label 1 is not used, we must change it to 1 and change every other labels
    if facets_set[0][0] != 1:
        a = facets_set[0][0]
        for i in range(len(facets_set)):
            for j in range(n):
                facets_set[i][j] -= (a - 1)
    m = max([facet[-1] for facet in facets_set])
    pic = m - n
    first_supposed_facet = list(range(1, n + 1))
    # if the first facet is not the regular n-simplex [1,...,n] then we relabel every faces to match it
    if facets_set[0] != first_supposed_facet:
        permutation = facets_set[0]
        for i in range(len(facets_set)):
            for j in range(n):
                a = facets_set[i][j]
                if a <= n:
                    facets_set[i][j] = permutation[a]
            facets_set[i].sort()
        facets_set.sort()
    # At this point, only the labels (n+1, ..., m) should be modified to obtain the minimal pure simplicial complex
    # for the lexicographic order
    minimal_facets_set = [facet.copy() for facet in facets_set]
    for last_labels_permutation_iter in permutations(range(n + 1, m + 1)):
        temp_facets_set = [facet.copy() for facet in facets_set]
        last_labels_permutation = list(last_labels_permutation_iter)
        for i in range(len(facets_set)):
            for j in range(n):
                a = temp_facets_set[i][j]
                if a > n:
                    temp_facets_set[i][j] = last_labels_permutation[a - n - 1]
            temp_facets_set[i].sort()
        temp_facets_set.sort()
        if temp_facets_set < minimal_facets_set:
            minimal_facets_set = [facet.copy() for facet in temp_facets_set]
    return minimal_facets_set



# setting up the connection with database
conn = sqlite3.connect("Real_Toric_Spaces.db")
print("Successfully connected to the database")

cursor = conn.cursor()

query = "INSERT OR IGNORE INTO PLSpheres (n, m, Pic, max_faces, min_non_faces, f_vector, h_vector, g_vector, " \
        "who_found_me, Date) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?); "

facets_list = read_file('./PLS_10_6')
# storing values in a variable
values = []

for k in range(len(facets_list)):
    print((k/len(facets_list)) * 100)
    facets = facets_list[k]
    facets_converted = json.loads(facets)
    n, m, pic = give_dim(facets_converted)
    faces_set = create_faces_set(facets_converted)
    f_vector = create_f_vector(faces_set)
    h_vector = create_h_vector(f_vector)
    g_vector = create_g_vector(h_vector)
    MNF_set = faces_set_to_MNF_set(faces_set)
    min_non_faces = []
    MNF_set.TreeToList(min_non_faces)
    author = 'Hyuntae Jang'
    date = datetime.datetime.now()
    values.append((n, m, pic, facets, str(min_non_faces), str(f_vector), str(h_vector), str(g_vector), author,
                   date.strftime("%x,%X")))
# executing the query with values

cursor.executemany(query, values)
# to make final output we have to run the 'commit()' method of the database object

conn.commit()
print(cursor.rowcount, "records inserted")
# closing the connection
cursor.close()
conn.close()
