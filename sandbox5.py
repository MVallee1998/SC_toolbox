import SimplicialComplex as sc
import json
import numpy as np
from itertools import combinations

def launch_all_lifts(n,k):
    m = n+k
    M = np.zeros((n, m))
    M[:, :n] = np.eye(n)
    list_vectors = []
    for j in range(1, 2**n):
        if j in sc.list_2_pow: continue
        vector = np.zeros(n)
        for i in range(n):
            if (j >> i) & 1:
                vector[i] = 1
        list_vectors.append(vector.copy())
    for iter in combinations(list_vectors,k):
        l = n
        for vector in iter:
            M[:,l] = vector.T
            l+=1
        if not sc.lifting_binary_matroid(M):
            print(M,"no lift")

launch_all_lifts(4,11)
