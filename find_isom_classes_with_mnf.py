import SimplicialComplex as sc
import json
m=10
n=6

intermediate_results_path = 'result/PLS_%d_%d_lin_alg_good_seeds' % (m, n)
final_results_path = 'final_results_BAK/PLS_%d_%d_2' % (m, n)

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


results = [json.loads(facets_bytes) for facets_bytes in read_file(intermediate_results_path)]

N = len(results)
K1 = sc.PureSimplicialComplex(results[0])
eq_classes = [K1]
for i in range(1,N):
    if i%1000 == 0:
        print((i/N)*1000,"%",len(eq_classes))
    K2 = sc.PureSimplicialComplex(results[i])
    is_isom = False
    for K1 in eq_classes:
        if sc.are_isom(K1,K2):
            is_isom = True
            break
    if not is_isom:
        eq_classes.append(K2)

data_to_text = []
for K in eq_classes:
    data_to_text.append(K.facets_bin)

text(data_to_text,final_results_path)