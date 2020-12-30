import binary_search_tree_global as bst
import SimplicialComplex as sc
import linear_alg_method as lam


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
    # Here we will try every face of dimension k
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

# char_funct = [6, 9, 10, 12, 7, 11, 13, 1
# 4, 15,8, 4, 2, 1]
# m = 13
# n = 9
# M, facets, ridges = lam.construct_matrix(char_funct, n, m)
# K = []
#
# for face in facets:
#     K.append(sc.binary_to_face(face,m))
# faces_set = create_faces_set(K)
# result = []
# faces_set_to_MNF_set(faces_set).TreeToList(result)
# print(result)