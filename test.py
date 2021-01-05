import MNF_set as mnf
import SimplicialComplex as sc
import json
import timeit
from itertools import combinations, permutations

m = 8
n = 4


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_good_seeds_lin_alg' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()

results = read_file('result/PLS_%d_%d_good_seeds_final' % (m, n))

for i in range(len(results)):
    K1_bin = json.loads(results[i])
    K1 = sc.PureSimplicialComplex(K1_bin)
    for j in range(i+1,len(results)):
        K2_bin = json.loads(results[j])
        K2 = sc.PureSimplicialComplex(K2_bin)
        print(sc.are_isom_to(K1,K2),K1.MNF_set,K2.MNF_set)