import SimplicialComplex as sc
import json
import timeit
import numpy as np
from itertools import combinations

m = 15
n = 11
p = m - n


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'final_results/PLS_%d_%d_sample' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


def decitoBin(numb):
    # checking if the given number is greater than 1
    if numb > 1:
        # if it is greater than 1 then use recursive approach by dividing number by 2
        decitoBin(numb // 2)
    # printing the binary representation of the given number
    print(numb % 2, end='')

results_path = 'final_results/PLS_%d_%d_sample'% (m, n)

for facets in [json.loads(facets_bytes) for facets_bytes in read_file(results_path)]:
    suspension_facets = [facet + sc.list_2_pow[14] for facet in facets]
    K = sc.PureSimplicialComplex(facets)
    K.compute_MNF_set()
    K.MNF_bin_to_MNF()
    print(K.MNF_set)


