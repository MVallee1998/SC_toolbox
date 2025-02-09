import SimplicialComplex as sc
import graph_methods as gm
from multiprocessing import Pool
from itertools import combinations, permutations
import timeit

n = 4
m = 7


def text(result):
    name = 'result/PLS_%d_%d_temp' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


def f(data):
    facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, starting_point = data
    result = gm.graph_method3_with_rec(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer,
                                       ridges_of_facets, n, m,
                                       1, -1,
                                       starting_point)
    return result


if __name__ == '__main__':
    final_result = []
    facets, ridges, G = gm.construct_complete_graph(n, m)
    facet_layer = gm.find_layer_v2(facets,ridges,m,0)
    facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets = gm.initialize_graph_method3(facets, ridges, G,
                                                                                                    facet_layer, n, m)

    step1 = gm.graph_method3_with_rec(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets,
                                      n, m, 0, 1)
    step1_good = []
    for starting_point in step1:
        step1_good.append(
            (facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, starting_point))
    print(len(step1_good))
    with Pool(processes=18) as pool:  # start 4 worker processes
        big_result = pool.map(f, step1_good)
        for result in big_result:
            for K in result:
                if K not in final_result:
                    final_result.append(K)
    text(final_result)

#
#
# def f(starting_point):
#     result = sc.graph_method2(facets, ridges, G, 0, starting_point)
#     text(result)
#
#
# if __name__ == '__main__':
#     with Pool(processes=19) as pool:  # start 4 worker processes
#         pool.map(f, step2)
