import SimplicialComplex as sc
import json
import numpy as np


m = 7
n = 3

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

final_results_path = 'final_results_BAK/PLS_%d_%d' % (m, n)

results = [json.loads(facets_bytes) for facets_bytes in read_file(final_results_path)]

for K_MF in results:
    K = sc.PureSimplicialComplex(K_MF)
    list_IDCM = sc.IDCM_Garrison_Scott(K)
    for IDCM in list_IDCM:
        print(sc.find_Z4_homology(K,IDCM[:K.n]))

# Pentagon = sc.PureSimplicialComplex([],[[1,3],[1,4],[2,4],[2,5],[3,5]],2)
# Pentagon.compute_facets_from_MNF_set()
# pendant3 = sc.multiple_wedge(Pentagon,[2,0,0,0,0])
# list_CM_bin = sc.Garrison_Scott(pendant3)
# list_CM_1 = []
# print(pendant3.MNF_set)
# n = pendant3.n
# m = pendant3.m
# for CM_bin in list_CM_bin:
#     CM = np.zeros((n,m))
#     for j in range(m):
#         for i in range(n):
#             if CM_bin[j] | sc.list_2_pow[i] == CM_bin[j]:
#                 CM[i,j] = 1
#     print(CM)
