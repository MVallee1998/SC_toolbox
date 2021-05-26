import SimplicialComplex as sc
import timeit
import graph_methods as gm

list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1] * 2)


char_funct = [3, 5, 6, 7, 1,2,4]
m = 7
n = 4

candidate_facets, ridges = gm.enumerate_facets_and_ridges(char_funct, n, m)



K0 = sc.PureSimplicialComplex([list_2_pow[n] - 1])

results = []
pile = [K0]

start = timeit.default_timer()
sc.Hyuntae_algo(pile, candidate_facets[1:], results, m)
stop = timeit.default_timer()
print(stop - start)
print(results)
