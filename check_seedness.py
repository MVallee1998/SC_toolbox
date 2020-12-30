import MNF_set as mnf
import SimplicialComplex as sc
import json
import timeit
import collections

m = 6
n = 2


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_lin_alg_mnf' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = read_file('result/PLS_%d_%d_temp' % (m, n))
l=0
counter = 0
start = timeit.default_timer()
for K_bytes in results:
    l += 1
    if l % 50 == 0:
        stop = timeit.default_timer()
        print("time spent :", stop - start)
        start = timeit.default_timer()
        print(l / len(results))
    K_bin = json.loads(K_bytes)
    K_sp = sc.PureSimplicialComplex(K_bin)
    K_sp.faces_set_to_MNF_set()
    print(K_sp.m, K_sp.n,K_sp.MNF_set_bin)

