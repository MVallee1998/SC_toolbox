import SimplicialComplex as sc
from multiprocessing import Pool
from itertools import combinations, permutations


# tests = []
# for combi_iter in combinations([3,5,6,9,10,12,7,11,13,14,15],6):
#     tests.append(list(combi_iter)+[8,4,2,1])
# for perm_iter in permutations([3, 5, 6, 7]):
#     tests.append(list(perm_iter) + [4, 2, 1])




char_funct = [ 12, 7, 11, 13, 14, 15, 8, 4, 2, 1]
n = 6
m = 10
facets, ridges, G = sc.construct_graph(char_funct,n,m)


facet_layer = [-1 for facet in facets]
facet_layer[0] = 0
sc.find_layer(G, facet_layer, [0], 1)
facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets = sc.initialize_graph_method3(facets, ridges,G,
                                                                                             facet_layer,n,m)
def f(starting_point):
    print(sc.graph_method3(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n, m, 1, -1,
                     starting_point))



if __name__ == '__main__':
    step1 = sc.graph_method3(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n, m, 0,
                             1)
    print(len(step1))
    with Pool(processes=8) as pool:  # start 4 worker processes
        pool.map(f, step1)


# name = 'result/PLS_test'
# t = open(name, mode='a', encoding='utf-8')
# def text(result):
#     for K in result:
#         t.write(str(K) + '\n')
#
#
# def f(starting_point):
#     result = sc.graph_method2(facets, ridges, G, 0, starting_point)
#     text(result)
#
#
if __name__ == '__main__':
    with Pool(processes=19) as pool:  # start 4 worker processes
        pool.map(f, step2)