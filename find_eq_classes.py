import SimplicialComplex as sc
import json
import timeit
import sqlite3
from multiprocessing import Pool

m = 11
n = 7

intermediate_results_path = 'intermediate_results/PLS_%d_%d' % (m, n)
final_results_path = 'final_results/PLS_%d_%d' % (m, n)


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




results = [json.loads(facets_bytes) for facets_bytes in read_file(intermediate_results_path)]
N = len(results)
print(N)
# setting up the connection with database
conn = sqlite3.connect("computed_PL.db")
print("Successfully connected to the database")
cursor = conn.cursor()
K0 = sc.PureSimplicialComplex(results[0])
K0_facets_mini = K0.find_minimal_lexico_order_db(cursor)
eq_classes = [K0_facets_mini]
counter = 0
start_sub = timeit.default_timer()
start = timeit.default_timer()

for i in range(1, N):
    if i % 10 == 0:
        stop_sub = timeit.default_timer()
        print("time spent for 10:", stop_sub - start_sub, " Percentage processed: ", i / N * 100,
              "%")
        start_sub = timeit.default_timer()
    K1 = sc.PureSimplicialComplex(results[i])
    K1_facets_mini = K1.find_minimal_lexico_order_db(cursor)
    cursor.execute('SELECT * FROM PLSpheres WHERE facets_bin="%s"' % json.dumps(K1_facets_mini))
    if cursor.fetchone() and K1_facets_mini not in eq_classes:
        eq_classes.append(K1_facets_mini)
        print(K1_facets_mini, counter)
        counter += 1
stop = timeit.default_timer()
print("Seeds selected.", " Time spent:", stop - start)
final_results = []
conn.close()
for K in eq_classes:
    final_results.append(K)
text(final_results, final_results_path)
