import MNF_set as mnf
import SimplicialComplex as sc
import json
import timeit
from itertools import combinations, permutations

m = 10
n = 6


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

K = sc.PureSimplicialComplex([list(range(1,n+1))])
K.list_unclosed_ridges()
print(K.unclosed_ridges)
pile = [K]
results = []
candidate_facets = []
for facet_iter in combinations(range(1,m+1),n):
    candidate_facets.append(sc.face_to_binary(list(facet_iter),m))
print(len(candidate_facets))
sc.Hyuntae_algo(pile,candidate_facets,results,m)
print(results)