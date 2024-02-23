import SimplicialComplex as sc
import json
import timeit
import tqdm

m = 10
n = 6
p = m - n
raw_results_path = 'raw_results/weak_psdmfd_%d_%d' % (m, n)
final_results_path = 'final_results/weak_psdmfd_%d_%d' % (m, n)


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

eq_classes = []
# for facets_isom in [json.loads(facets_bytes) for facets_bytes in read_file(initial_isom_path)]:
#     eq_classes.append(sc.PureSimplicialComplex(facets_isom))
for i in tqdm.tqdm(range(N)):
    K2 = sc.PureSimplicialComplex(results[i])
    # if not K2.is_a_seed():
    #     continue
    is_isom = False
    for K1 in eq_classes:
        if sc.are_isom(K1, K2):
            is_isom = True
            break
    if not is_isom and K2.n==n and K2.m==m:
        eq_classes.append(K2)
    else:
        del K2

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
