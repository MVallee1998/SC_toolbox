import SimplicialComplex as sc
import numpy as np
import timeit
import graph_methods as gm

# char_funct = [3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15,8,4,2,1]
# n = 11
# m = 15

# char_funct = [7, 11, 13, 14, 15, 8, 4, 2, 1]
# m = 9
# n = 5

char_funct = [3, 5, 6, 7, 4, 2, 1]
m = 7
n = 4

facets, ridges, G = gm.construct_graph(char_funct, n, m)
facet_layer = -1 * np.ones(len(facets))
facet_layer[0] = 0
gm.find_layer(G, facet_layer, [0], 1)
facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets = gm.initialize_graph_method3(facets, ridges, G,
                                                                                                facet_layer, n, m)
start = timeit.default_timer()
step1 = gm.graph_method3(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n, m, 0, -1)
print(len(step1))
stop = timeit.default_timer()
print("Time spent: ", stop - start)
