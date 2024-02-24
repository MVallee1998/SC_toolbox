import SimplicialComplex as sc
import json
import timeit
import tqdm

# raw_results_path = 'results_Hyuntae/PLS_10_6'
raw_results_path = 'results_Hyuntae/PLS_10_6'
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
# print(results[0])
# table = [0]*3
isom = []

N = len(results)

# i0 = 0
# K_i0 = sc.PureSimplicialComplex(results[i0])
# while not K_i0.is_a_seed() or not K_i0.is_promising() or not K_i0.is_Z2_homology_sphere() or not K_i0.is_closed():
#     i0 += 1
#     K_i0 = sc.PureSimplicialComplex(results[i0])
# K1 = sc.PureSimplicialComplex(results[i0])
eq_classes = []
# for facets_isom in [json.loads(facets_bytes) for facets_bytes in read_file(initial_isom_path)]:
#     eq_classes.append(sc.PureSimplicialComplex(facets_isom))
start_sub = timeit.default_timer()
start = start_sub
for i in tqdm.tqdm(range(N)):
    stop_sub = timeit.default_timer()
    K2 = sc.PureSimplicialComplex(results[i])
    K2.compute_MNF_set()
    # print(K2.MNF_set_bin)
    K2.MNF_bin_to_MNF()
    # print(K2.MNF_set)
    # if not K2.is_a_seed():
    #     continue
    is_isom = False
    # print("Time spent for 1", stop_sub - start_sub, (i / N) * 100, "%", len(eq_classes))
    start_sub = timeit.default_timer()
    a,b,c = K2.is_Z2_homology_sphere(), K2.is_promising(), K2.is_closed()
    if a and b and c and K2.m==10:
        eq_classes.append(K2)
    else:
        print(a,b,c)
        del K2
stop = timeit.default_timer()
print(len(eq_classes), " isomorphic classes found", " Time spent:", stop - start)
# for facets in results:
#     K = sc.PureSimplicialComplex(facets)

#     if K.m==8:
#         table[0]+=1
#     if K.m == 9:
#         table[1]+=1
#     if K.m == 10:
#         table[2]+=1
# print(table)
