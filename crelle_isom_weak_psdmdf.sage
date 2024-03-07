import json
from tqdm import tqdm
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


list_2_pow = [1]
for k in range(30):
    list_2_pow.append(list_2_pow[-1] * 2)

def face_to_binary(face, m):
    binary_face = m * [0]
    for k in face:
        binary_face[m - k] = 1
    return int(''.join(map(str, binary_face)), 2)

def binary_to_face(x, m):
    return [k + 1 for k in range(m) if (list_2_pow[k] | x) == x]

raw_results_path = 'results_crelle/seed_weak_psdmfd_%d_%d' % (m,n)
results = [json.loads(facets_bytes) for facets_bytes in read_file(raw_results_path)]
N=len(results)
print("imported")
isom_weak_psdmfd = []
for facets_bin in tqdm(results):
    facets = [binary_to_face(facet_bin,m) for facet_bin in facets_bin]
    K1 = SimplicialComplex(facets)
    # K1_MNF = SimplicialComplex(K1.minimal_nonfaces())
    is_isom = False
    for K2 in isom_weak_psdmfd:
        #K2_MNF = SimplicialComplex(K1.minimal_nonfaces())
        if K2.is_isomorphic(K1):
            is_isom=True
            break
    if not is_isom:
        isom_weak_psdmfd.append(K1)
    else:
        del K1
print("third row: ",len(isom_weak_psdmfd))
Crelle_result_path = 'results_Crelle/seed_weak_psdmfd_isom_%d_%d' % (m,n)
text([face_to_binary([list(facet) for facet in K.facets()],m) for K in isom_weak_psdmfd],Crelle_result_path)




