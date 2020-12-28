import SimplicialComplex as sc
import json
import timeit
import collections

m = 7
n = 3


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_1' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = read_file('result/PLS_%d_%d_temp1' % (m, n))

dictionary = dict()
for K_bin in results:
    dictionary[json.dumps(json.loads(K_bin))] = False
counter = 0
for K_str in dictionary:
    K_bin = json.loads(K_str)
    print(dictionary.get(K_str))
    if not dictionary.get(K_str):
        K_sc = sc.PureSimplicialComplex(K_bin)
        if K_sc.Pic == 4 and K_sc.is_closed():
            K_mini = K_sc.find_minimal_lexico_order(dictionary)
            if K_mini <= K_bin and K_mini not in results:
                if K_sc.is_Z2_homology_sphere() and K_sc.is_promising():
                    counter += 1
                    print(K_mini, counter)
                    results.append(K_mini)

text(results)
