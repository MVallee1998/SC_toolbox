import SimplicialComplex as sc
import json
import numpy as np

m = 9
n = 5
p = m - n
db_path = 'final_results/CSPLS_%d_%d' % (n, m)


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

results = [json.loads(facets_bytes) for facets_bytes in read_file(db_path)]

def find_all_lifts(original_facets,char_map):
    tree = np.zeros((n, m))
    tree[:, :n + 1] = char_map[:, :n + 1]
    for k in range(n + 1, m):
        tree[(tree[:, n:k] == 0).all(axis=1), k] = char_map[(tree[:, n:k] == 0).all(axis=1), k]
        rows_having_an_edge_previously = (tree[:, n:k] == 1).any(axis=1)
        rows_having_an_edge_now = (char_map[:, k] == 1)
        if np.logical_and(rows_having_an_edge_previously, rows_having_an_edge_now).any():
            index = np.where(np.logical_and(rows_having_an_edge_previously, rows_having_an_edge_now))[0][0]
            tree[index, k] = 1
    undecided_edges = char_map - tree
    relabelling = char_map - undecided_edges


for facets_bin in results:
    K = sc.PureSimplicialComplex(facets_bin)
    K.create_f_vector()
    print(K.f_vector)
    for char_map_bin in sc.Garrison_Scott(K):
        print(sc.char_funct_bin_to_numpy(char_map_bin,n))