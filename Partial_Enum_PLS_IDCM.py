from itertools import combinations, permutations
import numpy as np
import timeit


def Partial_Enum_PLS(K_p, m, n):
    start = timeit.default_timer()
    P = []
    Facets_sets = []
    Faces_in_Facets, Facets = ind_Facets(m, n)
    Ridges_with_Faces, Ridges, Faces = ind_Ridges(m, n)
    Facets_with_Ridges = []
    K = []
    ind_K = []
    Vertices = []
    Candidates = []
    Single_Ridges_with_Faces = []
    Closed_Faces = []
    Closed_Vertices = []
    Single_Ridges = []
    for i in range(len(Facets)):
        Closed_Vertices.append(set())
        Facets_sets.append(set(Facets[i]))
        K.append([])
        ind_K.append([])
        Vertices.append(set())
        Candidates.append([])
        Single_Ridges_with_Faces.append([])
        Closed_Faces.append([])
        Single_Ridges.append([])
        for j in range(len(Ridges_with_Faces)):
            Single_Ridges_with_Faces[-1].append([])
            Closed_Faces[-1].append(set())
            for k in range(len(Ridges_with_Faces[j])):
                Single_Ridges_with_Faces[-1][-1].append(set())
        Closed_Faces[-1].append(set())
    for i in range(len(Ridges)):
        Facets_with_Ridges.append(set())
        for j in range(len(Facets)):
            for k in Ridges[i]:
                if k not in Facets[j]:
                    break
            else:
                if Ridges[i][n - 2] == Facets[j][n - 2]:
                    Facets_with_Ridges[i].add(j)
    K[0] = Facets[0]
    ind_K[0] = 0
    for i in range(len(Single_Ridges_with_Faces[0])):
        for j in Faces_in_Facets[0][i]:
            Single_Ridges_with_Faces[0][i][j] = Ridges_with_Faces[i][j] & Faces_in_Facets[0][-1]
    Vertices[0] = Facets_sets[0]
    Single_Ridges[0] = Faces_in_Facets[0][-1]

    for i in range(1, len(K_p)):
        ind_K[i] = Facets.index(K_p[i])
        K[i] = K_p[i]
        Candidates[i] = set()
        Pre_Candi = set()
        for R in range(min(Single_Ridges[i - 1]) + 1):
            Pre_Candi.update(Facets_with_Ridges[R])
        for F in Pre_Candi:
            if F <= ind_K[i - 1] or F <= ind_K[i]:
                continue
            b = 1
            for dim in range(len(Closed_Faces[i - 1])):
                for CF in Closed_Faces[i - 1][dim]:
                    if CF in Faces_in_Facets[F][dim]:
                        b = 0
                        break
                if b == 0:
                    break
            else:
                Candidates[i].add(F)
        F = ind_K[i]
        Vertices[i] = Vertices[i - 1] | Facets_sets[F]
        for dim in range(len(Single_Ridges_with_Faces[i - 1])):
            for f in range(len(Single_Ridges_with_Faces[i - 1][dim])):
                Single_Ridges_with_Faces[i][dim][f] = Single_Ridges_with_Faces[i - 1][dim][f]
        for dim in range(1, n - 2):
            Closed_Faces[i][dim] = set(Closed_Faces[i - 1][dim])
            for f in Faces_in_Facets[F][dim]:
                Single_Ridges_with_Faces[i][dim][f] = Single_Ridges_with_Faces[i - 1][dim][f] ^ (
                            Ridges_with_Faces[dim][f] & Faces_in_Facets[F][-1])
                if Single_Ridges_with_Faces[i][dim][f] == set() and Single_Ridges_with_Faces[i - 1][dim][f] != set():
                    Closed_Faces[i][dim].add(f)
        intersection = Single_Ridges[i - 1] & Faces_in_Facets[F][-1]
        Closed_Faces[i][-1] = Closed_Faces[i - 1][-1] | intersection
        Single_Ridges[i] = (Single_Ridges[i - 1] | Faces_in_Facets[F][-1]) - intersection
        dim = 0
        Closed_Vertices[i] = set()
        Closed_Faces[i][dim] = set(Closed_Faces[i - 1][dim])
        for f in Faces_in_Facets[F][dim]:
            Single_Ridges_with_Faces[i][dim][f] = Single_Ridges_with_Faces[i - 1][dim][f] ^ (
                        Ridges_with_Faces[dim][f] & Faces_in_Facets[F][-1])

            if Single_Ridges_with_Faces[i][dim][f] == set() and Single_Ridges_with_Faces[i - 1][dim][f] != set():
                Closed_Faces[i][dim].add(f)
                Closed_Vertices[i].add(f + 1)
        Closed_Vertices[i] = Closed_Vertices[i] | Closed_Vertices[i - 1]
    if Single_Ridges[i] == set():
        c = 0
    else:
        c = 1
    i = i + 1

    # here is the main algorithm 4.4
    while i > 0:
        if c == 1:
            Candidates[i] = set()
            Pre_Candi = set()
            for R in range(min(Single_Ridges[i - 1]) + 1):
                Pre_Candi.update(Facets_with_Ridges[R])
            for F in Pre_Candi:
                if F <= ind_K[i - 1]:
                    continue
                b = 1
                for dim in range(len(Closed_Faces[i - 1])):
                    for CF in Closed_Faces[i - 1][dim]:
                        if CF in Faces_in_Facets[F][dim]:
                            b = 0
                            break
                    if b == 0:
                        break
                else:
                    Candidates[i].add(F)
        c = 1
        if len(Candidates[i]) == 0:
            i = i - 1
            c = 0
            continue
        F = min(Candidates[i])
        Candidates[i].remove(F)
        Vertices[i] = Vertices[i - 1] | Facets_sets[F]
        if max(Vertices[i]) != len(Vertices[i]):
            c = 0
            continue
        ind_K[i] = F
        K[i] = Facets[F]
        K_temp = K[:i + 1]
        New_CF = []
        for dim in range(len(Single_Ridges_with_Faces[i - 1])):
            for f in range(len(Single_Ridges_with_Faces[i - 1][dim])):
                Single_Ridges_with_Faces[i][dim][f] = Single_Ridges_with_Faces[i - 1][dim][f]
        for dim in range(1, n - 2):
            New_CF.append([])
            Closed_Faces[i][dim] = set(Closed_Faces[i - 1][dim])
            for f in Faces_in_Facets[F][dim]:
                Single_Ridges_with_Faces[i][dim][f] = Single_Ridges_with_Faces[i - 1][dim][f] ^ (
                            Ridges_with_Faces[dim][f] & Faces_in_Facets[F][-1])
                if Single_Ridges_with_Faces[i][dim][f] == set() and Single_Ridges_with_Faces[i - 1][dim][f] != set():
                    Closed_Faces[i][dim].add(f)
                    New_CF[-1].append(f)

        intersection = Single_Ridges[i - 1] & Faces_in_Facets[F][-1]
        Closed_Faces[i][-1] = Closed_Faces[i - 1][-1] | intersection
        Single_Ridges[i] = (Single_Ridges[i - 1] | Faces_in_Facets[F][-1]) - intersection
        dim = 0
        Closed_Vertices[i] = set()
        Closed_Faces[i][dim] = set(Closed_Faces[i - 1][dim])

        for f in Faces_in_Facets[F][dim]:
            Single_Ridges_with_Faces[i][dim][f] = Single_Ridges_with_Faces[i - 1][dim][f] ^ (
                        Ridges_with_Faces[dim][f] & Faces_in_Facets[F][-1])

            if Single_Ridges_with_Faces[i][dim][f] == set() and Single_Ridges_with_Faces[i - 1][dim][f] != set():
                Closed_Faces[i][dim].add(f)
                Closed_Vertices[i].add(f + 1)
        t = 1
        for d in range(1, n - 2):
            for f in New_CF[-d]:
                Link_face = K_temp
                new_cf = Faces[n - 2 - d][f]
                for l in range(len(new_cf)):
                    Link_face = Lk(Link_face, new_cf[l])
                Chain, Ind = Chain_cpx(Link_face, n - len(new_cf))
                if is_Homology_Sphere(Chain, Ind) == 0:
                    c = 0
                    t = 0
                    break
            else:
                continue
            break
        if t == 0:
            continue

        if is_Combi_for_Vert(K_temp, Closed_Vertices[i], n) == 0:
            c = 0
            continue
        Chain, Ind = Chain_cpx(K_temp, n)
        Closed_Vertices[i] = Closed_Vertices[i] | Closed_Vertices[i - 1]

        if is_lexico_min(K_temp, Vertices[i], Closed_Vertices[i], Chain, Ind, n) == 0:
            c = 0
            continue
        if is_IDCM(K_temp, m, n) == 0:
            c = 0
            continue

        if Single_Ridges[i] == set():
            c = 0
            if is_Homology_Sphere(Chain, Ind) == 1:
                P.append(K_temp)
                # text(K_temp,m,n)
        else:
            i += 1
    stop = timeit.default_timer()
    print(stop - start)
    return P


def is_IDCM(K, m, n): #is inductive characterisric map
    V = set(range(1, m + 1))
    C = []
    for i in K:
        C.append(list(V - set(i)))
    FP = Face_Poset(C)
    L = Z2n(m - n)
    M = [0] * n + L[0:m - n]
    L = L[m - n:]
    S = [0] * n
    i = n - 1
    c = 1
    while i < n:
        if c == 1:
            S[i] = L[:]
            for j in range(i, n - 1):
                S[i].remove(M[j + 1])
            for j in range(len(FP)):
                if i + 1 == FP[j][0]:
                    s = [0] * (m - n)
                    for k in range(1, len(FP[j])):
                        for l in range(len(s)):
                            if M[FP[j][k] - 1][l] == 1:
                                s[l] = (s[l] + 1) % 2
                    try:
                        S[i].remove(s)
                    except:
                        pass
        c = 1
        if S[i] == []:
            i = i + 1
            c = 0
            continue
        l = S[i][0]
        M[i] = l
        S[i].remove(l)
        if i == 0:
            return 1
        else:
            i += -1
    return 0


def Z2n(n):
    V = []
    N = range(n)
    for i in range(1, n + 1):
        C = combinations(N, i)
        for c in C:
            v = []
            for j in range(n):
                if j in c:
                    v.append(1)
                    continue
                v.append(0)
            V.append(v)
    return V


def Face_Poset(K):
    FP = []
    for k in K:
        k.sort()
        for i in range(2, len(k) + 1):
            C = combinations(k, i)
            for c in C:
                if list(c) not in FP:
                    FP.append(list(c))
    return FP


def ind_Ridges(m, n):
    L = []
    Ridges = [list(c) for c in combinations(range(1, m + 1), n - 1)]
    Faces = []
    for i in range(1, n - 1):
        L.append([])
        Faces.append([])
        C = combinations(range(1, m + 1), i)
        for c in C:
            L[-1].append(set())
            Faces[-1].append(list(c))
            for R in range(len(Ridges)):
                for j in c:
                    if j not in Ridges[R]:
                        break
                else:
                    L[-1][-1].add(R)
    return L, Ridges, Faces


def ind_Facets(m, n):
    L = []
    Facets = []
    C = combinations(range(1, m + 1), n)
    for c in C:
        Facets.append(list(c))
        L.append([])
        for i in range(1, n):
            L[-1].append(set())
    for i in range(1, n):
        C = combinations(range(1, m + 1), i)
        ind = 0
        for c in C:
            for F in range(len(Facets)):
                for j in c:
                    if j not in Facets[F]:
                        break
                else:
                    L[F][i - 1].add(ind)
            L[-1][-1] = set(L[-1][-1])
            ind += 1
    return L, Facets


def Chain_cpx(K, n):  # ordered K
    C = [K[:]]
    ind = [[]]
    i = n - 1
    while i > 0:
        C.append([])
        ind.append([])
        j = n - i
        for k in range(len(C[j - 1])):
            D = combinations(C[j - 1][k], i)
            for d in D:
                c = list(d)
                if c in C[j]:
                    ind[j][C[j].index(c)].append(k)
                else:
                    C[j].append(c)
                    ind[j].append([k])
        i = i - 1
    return C, ind


def Z2_rank(M):
    m = len(M)
    n = len(M[0])
    i = 0
    j = 0
    while i < m and j < n:
        c = 0
        for p in range(i, m):
            if M[p, j] == 1:
                c = 1
                M[[i, p], :] = M[[p, i], :]
                for k in range(p + 1, m):
                    if M[k, j] == 1:
                        M[k, :] = (M[k, :] + M[i, :]) % 2
                break
        if c == 0:
            j += 1
            continue
        i += 1
        j += 1
    return i


def Z2_Homology(Chain, Ind):
    H = []
    im = 0
    for i in range(1, len(Chain)):
        M = np.zeros((len(Chain[i]), len(Chain[i - 1])))
        for j in range(len(Ind[i])):
            for k in range(len(Ind[i][j])):
                M[j, Ind[i][j][k]] = 1
        ker = len(Chain[i - 1]) - Z2_rank(M)
        H.append(ker - im)
        im = len(Chain[i - 1]) - ker
    H.append(len(Chain[-1]) - im)
    return H


def is_Homology_Sphere(Chain, Ind):
    H = Z2_Homology(Chain, Ind)
    if H == [2]:
        return 1
    L = [1] + [0] * (len(Chain) - 2) + [1]
    if H == L:
        return 1
    else:
        return 0


def is_Combi_for_Vert(K, V, n):
    for v in V:
        Chain, Ind = Chain_cpx(Lk(K, v), n - 1)
        if is_Homology_Sphere(Chain, Ind) == 0:
            return 0
    return 1


def Lk(K, v):
    L = []
    for i in K:
        if v not in i:
            continue
        L.append([])
        for j in i:
            if j != v:
                L[-1].append(j)
    return L


def relabel(K, V, NV):  # V, NV: lists
    NK = []
    for i in range(len(K)):
        NK.append([])
        for j in range(len(K[i])):
            NK[i].append(NV[V.index(K[i][j])])
        NK[i].sort()
    NK.sort()
    return NK


def is_lexico_min(K, V_set, CV, Chain, Ind, n):
    dCV = []
    d = len(Ind[-1][(Chain[-1].index([1]))])
    for v in CV:
        deg = len(Ind[-1][(Chain[-1].index([v]))])
        if deg < d:
            return 0
        elif deg == d:
            dCV.append(v)
    V = list(V_set)
    V.sort()
    NV = [0] * len(Chain[-1])
    for v in dCV:
        NV[0] = v
        for i in range(len(Chain[1])):
            if v in Chain[1][i]:
                p_set = set(Chain[1][i]) - {v}
                q_set = set(Chain[0][Ind[1][i][0]]) ^ set(Chain[0][Ind[1][i][1]])
                P = permutations(p_set)
                for p in P:
                    for j in range(len(p)):
                        NV[1 + j] = p[j]
                    Q = permutations(q_set)
                    for q in Q:
                        for k in range(len(q)):
                            NV[n - 1 + k] = q[k]
                        R = permutations(V_set - {v} - p_set - q_set)
                        for r in R:
                            for k in range(len(r)):
                                NV[n + 1 + k] = r[k]
                            if relabel(K, NV, V) < K:
                                return 0
    return 1


def text(K, m, n):
    name = 'PLS_%d_%d' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    t.write(str(K) + '\n')
    t.close()


def Fix_Pic(P, pic):
    L = []
    n = len(P[0][0])
    m = n + pic
    for i in P:
        M = 0
        for j in i:
            M = max(max(j), M)
        if M == m:
            L.append(i)
    return L

K_bin = [2047, 6111, 9727, 10223, 10235, 13791, 14287, 14299, 14329, 14330, 15865, 15866, 16121, 16122, 16249, 16250, 17407, 18175, 18303, 18367, 18423, 18429, 18430, 21471, 22239, 22367, 22431, 22487, 22493, 22494, 22521, 22522, 24057, 24058, 24313, 24314, 24441, 24442, 25087, 25583, 25595, 25855, 25983, 26047, 26103, 26109, 26110, 26351, 26363, 26479, 26491, 26543, 26555, 26599, 26605, 26606, 26611, 28153, 28154, 28409, 28410, 28537, 28538, 29151, 29647, 29659, 29689, 29690, 29919, 30047, 30111, 30167, 30173, 30174, 30415, 30427, 30543, 30555, 30607, 30619, 30649, 30650, 30663, 30669, 30670, 30675, 30697, 30698, 30705, 30706, 31225, 31226, 31481, 31482, 31609, 31610, 32185, 32186, 32217, 32218, 32233, 32234, 32241, 32242, 32441, 32442, 32473, 32474, 32489, 32490, 32497, 32498, 32569, 32570, 32601, 32602, 32617, 32618, 32625, 32626]
list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1]*2)
K= []
for facet_bin in K_bin:
    K.append([])
    for i in range(len(list_2_pow)):
        element = list_2_pow[i]
        if element | facet_bin == facet_bin:
            K[-1].append(i+1)
K.sort()
Chain, Ind = Chain_cpx(K,11)
print(is_Homology_Sphere(Chain,Ind))