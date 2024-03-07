import SimplicialComplex as sc
import json
import timeit
import tqdm

m=11
n=7

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
number_of_m = [0]*(m+1)
raw_results_path = 'results_crelle/seed_weak_psdmfd_%d_%d' % (m,n)
results = [json.loads(facets_bytes) for facets_bytes in read_file(raw_results_path)]
Crelle_result_path = 'results_Crelle/seed_weak_psdmfd_%d_%d' % (m,n)
isom_weak_psdmfd = []
for facets in tqdm.tqdm(results):
    K1 =sc.PureSimplicialComplex(facets)
    is_isom = False
    for K2 in isom_weak_psdmfd:
        if sc.are_isom(K1,K2):
            is_isom=True
            break
    if not is_isom:
        isom_weak_psdmfd.append(K1)
print("third row: ",len(seed_PLS))
Crelle_result_path = 'results_Crelle/seed_weak_psdmfd_isom_%d_%d' % (m,n)
text([list(K.facets_bin) for K in isom_weak_psdmfd],Crelle_result_path)




