import SimplicialComplex as sc
import graph_methods as gm
from multiprocessing import Pool
from itertools import combinations, permutations
import timeit

n = 6
m = 10
tests = []

def text(result):
    name = 'result/PLS_%d_%d_temp' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


def f(data):
    print("ninja")
    start = timeit.default_timer()
    facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, starting_point = data
    result = gm.graph_method3_with_rec(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n, m,
                              1, -1,
                              starting_point)
    stop = timeit.default_timer()
    if stop-start> 50:
        print(stop-start)
    return result


if __name__ == '__main__':
    k = 0
    final_result = []
    sizes = []
    start = timeit.default_timer()
    list_char_funct = sc.enumerate_char_funct_orbits(n, m)
    for char_funct in sorted(list_char_funct):
        part_result = []
        start_sub = timeit.default_timer()
        k += 1
        # for perm_iter in permutations([3, 5, 6, 7]):
        #     tests.append(list(perm_iter) + [4, 2, 1])
        print((k/len(list_char_funct))*100, char_funct)
        facets, ridges, G = gm.construct_graph(char_funct, n, m)
        facet_layer = [-1 for facet in facets]
        facet_layer[0] = 0
        gm.find_layer(G, facet_layer, [0], 1)
        facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets = gm.initialize_graph_method3(facets, ridges,
                                                                                                        G,
                                                                                                        facet_layer, n,
                                                                                                        m)

        step1 = gm.graph_method3_with_rec(facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, n,
                                 m, 0, 1)
        step1_good = []
        for starting_point in step1:
            step1_good.append(
                (facets, ridges, facets_for_ridges, facets_for_ridges_with_layer, ridges_of_facets, starting_point))
        print(len(step1_good))
        with Pool(processes=18) as pool:  # start 4 worker processes
            big_result = pool.map(f, step1_good)
            for result in big_result:
                for K in result:
                    # K_sp = sc.PureSimplicialComplex(K)
                    # if K_sp.Pic == 4 and K_sp.is_promising() and K_sp.is_Z2_homology_sphere():
                    #     part_result.append(K_sp)
                    if K not in final_result:
                        final_result.append(K)
            if len(part_result) not in sizes:
                sizes.append(len(part_result))
        stop_sub = timeit.default_timer()
        print("Time spent: ", stop_sub - start_sub)
    stop = timeit.default_timer()
    print("Time spent: ", stop - start)
    text(final_result)
    print(sizes)
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
