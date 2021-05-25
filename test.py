import SimplicialComplex as sc
import json

m = 12
n = 8

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

partial_results_path = 'partial_results/PLS_%d_%d' % (m, n)

results = [json.loads(facets_bytes) for facets_bytes in read_file(partial_results_path)]

for K_MF in results:
    K = sc.PureSimplicialComplex(K_MF)
    sc.IDCM_Garrison_Scott(K)