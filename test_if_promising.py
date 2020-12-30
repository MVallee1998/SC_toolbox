import SimplicialComplex as sc
import json
import timeit


m = 9
n = 5


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_lin_alg_test' % (m,n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = read_file('result/PLS_%d_%d_lin_alg' % (m, n))

K_result = []
l = 0
start = timeit.default_timer()
counter = 0
for K_bytes in results:
    l += 1
    if l % 50 == 0:
        stop = timeit.default_timer()
        print("time spent :", stop - start)
        start = timeit.default_timer()
        print(l / len(results))
    K_bin = json.loads(K_bytes)
    K_sp = sc.PureSimplicialComplex(K_bin)
    if K_sp.Pic == 4 and K_sp.is_promising() and K_sp.is_Z2_homology_sphere():
        K_sp_mini = K_sp.find_minimal_lexico_order()
        if K_sp_mini <= K_sp.facets_bin:
            K_bin = K_sp_mini
            if K_bin not in K_result:
                print(K_bin, counter)
                counter+=1
                K_result.append(K_bin)

    # list_2_pow = [1]
    # for k in range(8):
    #     list_2_pow.append(list_2_pow[-1]*2)
    # K= []
    # for facet_bin in K_bin:
    #     K.append([])
    #     for i in range(len(list_2_pow)):
    #         element = list_2_pow[i]
    #         if element | facet_bin == facet_bin:
    #             K[-1].append(i+1)
    # K.sort()
    # Chain, Ind = Partial_Enum_PLS_IDCM.Chain_cpx(K, 6)
    # if Partial_Enum_PLS_IDCM.is_Homology_Sphere(Chain, Ind) and K not in K_result:
    #     K_result.append(K)
K_result.sort()
text(K_result)
