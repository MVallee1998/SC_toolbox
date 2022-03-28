import SimplicialComplex as sc
import json
import timeit

m = 15
n = 11
p = m - n
raw_results_path = 'partial_results/PLS_11_15_seeds_test'
# initial_isom_path = 'final_results/PLS_14_10_1'
final_results_path = 'final_results/PLS_%d_%d_sample' % (m, n)


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result, path):
    t = open(path, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()


results = [json.loads(facets_bytes) for facets_bytes in read_file(raw_results_path)]

N = len(results)

i0 = 0
K_i0 = sc.PureSimplicialComplex(results[i0])
while not K_i0.is_promising() or not K_i0.is_Z2_homology_sphere() or not K_i0.is_closed():
    i0 += 1
    K_i0 = sc.PureSimplicialComplex(results[i0])
K1 = sc.PureSimplicialComplex(results[i0])
eq_classes = [K1]
K1.compute_MNF_set()
K1.MNF_bin_to_MNF()
print(K1.MNF_set)
# for facets_isom in [json.loads(facets_bytes) for facets_bytes in read_file(initial_isom_path)]:
#     eq_classes.append(sc.PureSimplicialComplex(facets_isom))
start_sub = timeit.default_timer()
start = start_sub
for i in range(N):
    stop_sub = timeit.default_timer()
    K2 = sc.PureSimplicialComplex(results[i])
    is_isom = False
    if (not K2.Pic == p) or (not K2.is_a_seed()):
        del K2
        continue
    print("Time spent for 1", stop_sub - start_sub, (i / N) * 100, "%", len(eq_classes))
    start_sub = timeit.default_timer()
    for K1 in eq_classes:
        if sc.are_isom(K1, K2):
            is_isom = True
            break
    if not is_isom and K2.is_Z2_homology_sphere() and K2.is_promising() and K2.is_closed():
        eq_classes.append(K2)
        K2.compute_MNF_set()
        K2.MNF_bin_to_MNF()
        print(K2.MNF_set)
    else:
        del K2
stop = timeit.default_timer()
print(len(eq_classes), " isomorphic classes found", " Time spent:", stop - start)

# N = len(eq_classes)
#
# good_seeds = []
# start = timeit.default_timer()
# start_sub = timeit.default_timer()
# for i in range(N):
#     K = eq_classes[i]
#     if i % 100 == 0:
#         stop_sub = timeit.default_timer()
#         print("time spent for 100:", stop_sub - start_sub, " Percentage processed: ", i / N * 100, "%")
#         start_sub = timeit.default_timer()
#     if K.is_promising() and K.is_closed() and K.is_Z2_homology_sphere():
#         good_seeds.append(K)
# stop = timeit.default_timer()
# print(len(good_seeds), "Good seeds selected", " Time spent:", stop - start)

data_to_text = []
for K in eq_classes:
    data_to_text.append(K.facets_bin)

text(data_to_text, final_results_path)
stop = timeit.default_timer()
print("Final result saved.", " Time spent:", stop - start)
