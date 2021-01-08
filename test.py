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

results = read_file('result/PLS_%d_%d_good_seeds_lin_alg_final' % (m, n))

for K_bin in [[[1, 2, 3, 4, 5], [1, 2, 3, 4, 6], [1, 2, 3, 5, 6], [1, 2, 4, 5, 6], [1, 3, 4, 5, 6], [2, 3, 4, 5, 7], [2, 3, 4, 6, 7], [2, 3, 5, 6, 7], [2, 4, 5, 6, 8], [2, 4, 5, 7, 8], [2, 4, 6, 7, 8], [2, 5, 6, 7, 8], [3, 4, 5, 6, 9], [3, 4, 5, 7, 8], [3, 4, 5, 8, 9], [3, 4, 6, 7, 9], [3, 4, 7, 8, 9], [3, 5, 6, 7, 9], [3, 5, 7, 8, 9], [4, 5, 6, 8, 9], [4, 6, 7, 8, 9], [5, 6, 7, 8, 9]]]:
    # K_bin = json.loads(K_byte)
    K_sc = sc.PureSimplicialComplex(K_bin)
    K_sc.compute_MNF_set()
    K_sc.MNF_bin_to_MNF()
    print(K_sc.MNF_set)
    K_sc.create_FP()
    for vertex in K_sc.FP_bin[0]:
        link_of_vertex = sc.Link_of(K_sc,vertex)
        link_of_vertex_mini_facets = link_of_vertex.find_minimal_lexico_order()
        link_of_vertex_mini_sc = sc.PureSimplicialComplex(link_of_vertex_mini_facets)
        link_of_vertex_mini_sc.compute_MNF_set()
        link_of_vertex_mini_sc.MNF_bin_to_MNF()
        print("Link of",sc.binary_to_face(vertex,K_sc.m), link_of_vertex_mini_sc.Pic, link_of_vertex_mini_sc.is_a_seed(), link_of_vertex_mini_sc.MNF_set)