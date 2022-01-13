import SimplicialComplex as sc
import json
import timeit
import collections

m = 9
n = 5


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_lin_alg_good_seeds_final' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = read_file('result/PLS_%d_%d_lin_alg_good_seeds' % (m, n))

dictionary = dict()
for K_bin in results:
    dictionary[json.dumps(json.loads(K_bin))] = False
counter = 0
l=0
taille = len(dictionary)
PL_Spheres = []
start =  timeit.default_timer()
for K_str in dictionary:
    l+=1
    if l % 100 == 0:
        stop = timeit.default_timer()
        print("time spent :", stop - start)
        start = timeit.default_timer()
        print(l / taille)
    K_bin = json.loads(K_str)
    # print(K_bin)
    if not dictionary.get(K_str):
        K_sc = sc.PureSimplicialComplex(K_bin)
        if K_sc.Pic == 4 and K_sc.is_Z2_homology_sphere() and K_sc.is_promising():
            K_mini = K_sc.find_minimal_lexico_order(dictionary)
            if K_mini <= K_bin and K_mini not in PL_Spheres:
                counter += 1
                print(K_mini, counter)
                PL_Spheres.append(K_mini)

text(PL_Spheres)
