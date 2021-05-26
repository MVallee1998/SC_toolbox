import SimplicialComplex as sc
import json
import timeit

m = 6
n = 2


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_lin_alg_all_seeds' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = read_file('result/PLS_%d_%d_lin_alg' % (m, n))
l=0
counter = 0
start = timeit.default_timer()
print(len(results))
list_seeds = []
for K_bytes in results:
    l += 1
    if l % 1000 == 0:
        stop = timeit.default_timer()
        print("time spent :", stop - start)
        start = timeit.default_timer()
        print(l / len(results))
    K = json.loads(K_bytes)
    K_sc = sc.PureSimplicialComplex(K)
    if K_sc.is_a_seed():
        list_seeds.append(K)

text(list_seeds)
