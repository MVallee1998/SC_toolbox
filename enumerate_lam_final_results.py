import SimplicialComplex as sc
import json
import timeit
m = 11
n = 7

raw_results_path = 'raw_results/PLS_%d_%d' % (m, n)
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
list_of_seeds = []
counter = 0
start_sub = timeit.default_timer()
start = timeit.default_timer()
for i in range(len(results)):
    counter+=1
    if counter % 1000 == 0:
        stop_sub = timeit.default_timer()
        print("time spent for 1000:", stop_sub - start_sub, " Percentage processed: ", counter / len(results) * 100, "%")
        start_sub = timeit.default_timer()
    facets = results[i]
    K = sc.PureSimplicialComplex(facets)
    if K.Pic == 4 and K.is_a_seed() and K not in list_of_seeds:
        list_of_seeds.append(facets)
    del K
stop = timeit.default_timer()
print("Seeds selected.", " Time spent:", stop - start)

list_of_good_seeds = []
counter = 0
start = timeit.default_timer()
start_sub = timeit.default_timer()
for i in range(len(list_of_seeds)):
    counter+=1
    if counter % 100 == 0:
        stop_sub = timeit.default_timer()
        print("time spent for 100:", stop_sub - start_sub, " Percentage processed: ", counter / len(list_of_seeds) * 100, "%")
        start_sub = timeit.default_timer()
    facets = list_of_seeds[i]
    K = sc.PureSimplicialComplex(facets)
    if K.Pic == 4 and K.is_promising() and K.is_Z2_homology_sphere() and K.is_closed():
        list_of_good_seeds.append(facets)
    del K
stop = timeit.default_timer()
print("Good seeds selected", " Time spent:", stop - start)

N = len(list_of_good_seeds)
list_final = []
counter = 0
counter_PLS = 0
dictionary_ref = dict()
start = timeit.default_timer()
start_sub = timeit.default_timer()

K1 = sc.PureSimplicialComplex(list_of_good_seeds[0])
eq_classes = [K1]
for i in range(1,N):
    if i%100 == 0:
        stop_sub = timeit.default_timer()
        print("Time spent for 100", stop_sub-start_sub,(i/N)*100,"%",len(eq_classes))
        start_sub = timeit.default_timer()
    K2 = sc.PureSimplicialComplex(list_of_good_seeds[i])
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
stop = timeit.default_timer()
print("Final result saved.", " Time spent:", stop - start)


