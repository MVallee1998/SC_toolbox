import SimplicialComplex as sc
import json
import timeit
import tqdm

m=10
n=6

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
raw_results_path = 'raw_results/weak_psdmfd_%d_%d' % (m,n)
results = [json.loads(facets_bytes) for facets_bytes in read_file(raw_results_path)]
chosed_m=[]
for facets in tqdm.tqdm(results):
    K =sc.PureSimplicialComplex(facets)
    number_of_m[K.m]+=1
    if K.m==m:
        chosed_m.append(K.facets_bin.copy())
    del K
print(list(range(m+1)))
del results
Crelle_result_path = 'results_Crelle/chosed_m_%d_%d' % (m,n)
text(chosed_m,Crelle_result_path)
# print(number_of_m)
print("first row: ",len(chosed_m))
seed_weak_psdmfd = []
for facets in tqdm.tqdm(chosed_m):
    K =sc.PureSimplicialComplex(facets)
    if K.is_a_seed():
        seed_weak_psdmfd.append(K.facets_bin.copy())
    del K
print("second row: ",len(seed_weak_psdmfd))
del chosed_m
Crelle_result_path = 'results_Crelle/seed_weak_psdmfd_%d_%d' % (m,n)
text(seed_weak_psdmfd,Crelle_result_path)
seed_PLS = []
for facets in tqdm.tqdm(seed_weak_psdmfd):
    K =sc.PureSimplicialComplex(facets)
    if K.is_Z2_homology_sphere() and K.is_closed() and K.is_promising():
        seed_PLS.append(K.facets_bin.copy())
    del K
print("third row: ",len(seed_PLS))
Crelle_result_path = 'results_Crelle/seed_PLS_%d_%d' % (m,n)
text(seed_PLS,Crelle_result_path)




