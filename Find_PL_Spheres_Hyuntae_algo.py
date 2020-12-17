import SimplicialComplex as sc
from multiprocessing import Pool
from itertools import combinations, permutations
import timeit

list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1] * 2)


char_funct = [7,11, 13, 14, 15, 8, 4, 2, 1]
m = 9
n = 5



K0 = sc.PureSimplicialComplex([list_2_pow[n] - 1])
ref_facets = []
for facet_iter in combinations(range(1, m + 1), n):
    ref_facets.append(sc.face_to_binary(list(facet_iter), m))
results = []
pile = [K0]

start = timeit.default_timer()
sc.Hyuntae_algo(pile, ref_facets[1:], results, m)
stop = timeit.default_timer()
print(stop - start)
print(results)
