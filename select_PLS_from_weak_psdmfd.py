import SimplicialComplex as sc
import json
import timeit
import tqdm
m = 9
n = 5
p = m-n
raw_results_path = 'final_results/weak_psdmfd_%d_%d' % (m, n)
final_results_path = 'final_results/PLS_%d_%d' % (m, n)

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
list_of_PLS = []


for i in tqdm.tqdm(range(len(results))):
    facets = results[i]
    K = sc.PureSimplicialComplex(facets)
    # K.compute_MNF_set()
    # K.MNF_bin_to_MNF()
    # print(K.MNF_set)
    if K.is_closed() and K.is_Z2_homology_sphere() and K.is_promising():
        list_of_PLS.append(K)
    else: del K

data_to_text = []
for K in list_of_PLS:
    data_to_text.append(K.facets_bin)

text(data_to_text,final_results_path)


