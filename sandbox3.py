import SimplicialComplex as sc
import json
import timeit
from itertools import combinations
import numpy as np

m = 15
n = 11
p = m - n



def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data



for n in range(2,12):
    m=n+4
    converted_path = 'converted/PLS_%d_%d' % (m, n)
    final_results_path = 'final_results/CSPLS_%d_%d' % (n, m)
    results = [json.loads(facets_bytes) for facets_bytes in read_file(final_results_path)]
    N = len(results)
    print(N)
    t = open(converted_path,'w')
    for k in range(N):
        K_facets_bin = results[k]
        K = sc.PureSimplicialComplex(K_facets_bin)
        t.write(str((np.array(K.facets)-1).tolist()).replace('[','{').replace(']','}').replace(' ','')+'\n')
    t.close()
