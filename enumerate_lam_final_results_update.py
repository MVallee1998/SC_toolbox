import SimplicialComplex as sc
import json
import timeit
m = 9
n = 4
p=m-n
raw_results_path = 'raw_results/PLS_%d_%d' % (m, n)
final_results_path = 'final_results/PLS_%d_%d_new3' % (m, n)

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

N = len(results)

i0=0
while (not sc.PureSimplicialComplex(results[i0]).is_a_seed()) or sc.PureSimplicialComplex(results[i0]).Pic != p:
    i0+=1
K1 = sc.PureSimplicialComplex(results[i0])
eq_classes = [K1]
start_sub = timeit.default_timer()
start = start_sub
for i in range(i0+1,N):
    if i%100 == 0:
        stop_sub = timeit.default_timer()
        print("Time spent for 100", stop_sub-start_sub,(i/N)*100,"%",len(eq_classes))
        start_sub = timeit.default_timer()
    K2 = sc.PureSimplicialComplex(results[i])
    is_isom = False
    if K2.Pic == p and K2.is_a_seed():
        for K1 in eq_classes:
            if sc.are_isom(K1,K2):
                is_isom = True
                break
        if not is_isom:
            eq_classes.append(K2)
        else: del K2
    else: del K2
stop = timeit.default_timer()
print(len(eq_classes)," isomorphic classes found", " Time spent:", stop - start)

N = len(eq_classes)

good_seeds = []
start = timeit.default_timer()
start_sub = timeit.default_timer()
for i in range(N):
    K = eq_classes[i]
    if i % 100 == 0:
        stop_sub = timeit.default_timer()
        print("time spent for 100:", stop_sub - start_sub, " Percentage processed: ", i / N * 100, "%")
        start_sub = timeit.default_timer()
    if K.is_promising() and K.is_closed() and K.is_Z2_homology_sphere():
        good_seeds.append(K)
stop = timeit.default_timer()
print(len(good_seeds),"Good seeds selected", " Time spent:", stop - start)

data_to_text = []
for K in good_seeds:
    data_to_text.append(K.facets_bin)

text(data_to_text,final_results_path)
stop = timeit.default_timer()
print("Final result saved.", " Time spent:", stop - start)


