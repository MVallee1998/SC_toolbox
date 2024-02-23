import SimplicialComplex as sc
import json
import timeit
m = 10
n = 6

raw_results_path = 'raw_results/weak_psdmfd_%d_%d' % (m, n)
final_results_path = 'final_results/PLS_%d_%d_test2' % (m, n)

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result,path):
    t = open(path, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = [json.loads(facets_bytes) for facets_bytes in read_file(raw_results_path)]
list_of_seeds_orbit = []
list_of_good_PLS = []
start_sub = timeit.default_timer()
start = timeit.default_timer()

for i in range(len(results)):
    if i % 1000 == 0:
        stop_sub = timeit.default_timer()
        print("time spent for 1000:", stop_sub - start_sub, " Percentage processed: ", i / len(results) * 100, "%")
        start_sub = timeit.default_timer()
    facets = results[i]
    K2 = sc.PureSimplicialComplex(facets)
    if K2.Pic == 4 and K2.is_a_seed():
        is_isom = False
        for K1 in list_of_seeds_orbit:
            start_isom = timeit.default_timer()
            if sc.are_isom(K1, K2):
                is_isom = True
                stop_isom = timeit.default_timer()
                if stop_isom - start_isom>2:
                    print("long MNF",K1.MNF_set)
                break
        if not is_isom:
            list_of_seeds_orbit.append(K2)
            if K2.is_promising() and K2.is_Z2_homology_sphere() and K2.is_closed():
                list_of_good_PLS.append(K2)
                K2.MNF_bin_to_MNF()
                print("New Good seed MNF",K2.MNF_set)
                print("Current number of PLS", len(list_of_good_PLS))


data_to_text = []
for K in list_of_good_PLS:
    data_to_text.append(K.facets_bin)

text(data_to_text,final_results_path)
stop = timeit.default_timer()
print("Final result saved.", " Time spent:", stop - start)


