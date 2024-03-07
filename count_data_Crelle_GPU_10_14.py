import SimplicialComplex as sc
import json
import timeit
import tqdm

m=14
n=10
k_max = 3


# raw_results_path = 'results_Hyuntae/PLS_10_6'

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
chosed_m=[]
for k in range(k_max):
    raw_results_path = 'raw_results/GPU_%d_%d_%d.out' % (n,m,k)
    results = [json.loads(facets_bytes) for facets_bytes in read_file(raw_results_path)]
    for facets in tqdm.tqdm(results):
        K =sc.PureSimplicialComplex(facets)
        number_of_m[K.m]+=1
        if K.m==m:
            chosed_m.append(K.facets_bin)
        del K
print(list(range(m+1)))
print(number_of_m)
del results
Crelle_result_path = 'results_Crelle/chosed_m_%d_%d' % (m,n)
text(chosed_m,Crelle_result_path)
# print(number_of_m)
print("first row: ",len(chosed_m))
# seed_weak_psdmfd = []
# for facets in tqdm.tqdm(chosed_m):
#     K =sc.PureSimplicialComplex(facets)
#     if K.is_a_seed():
#         seed_weak_psdmfd.append(K.facets_bin.copy())
#     del K
# print("second row: ",len(seed_weak_psdmfd))
# del chosed_m
# Crelle_result_path = 'results_Crelle/seed_weak_psdmfd_%d_%d' % (m,n)
# text(seed_weak_psdmfd,Crelle_result_path)
IDCM_PLS = []
for facets in tqdm.tqdm(chosed_m):
    K =sc.PureSimplicialComplex(facets)
    if K.is_promising() and K.is_Z2_homology_sphere() and K.is_closed():
        IDCM_PLS.append(list(K.facets_bin.copy()))
    del K
print("third row: ",len(IDCM_PLS))
Crelle_result_path = 'results_Crelle/IDCM_PLS_%d_%d' % (m,n)
text(IDCM_PLS,Crelle_result_path)






