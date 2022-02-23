import SimplicialComplex as sc
import json
import timeit
import numpy as np
from itertools import combinations

m = 15
n = 11
p = m-n
def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'partial_results/PLS_%d_%d_seeds' % (m, n)
    t = open(name, mode='a', encoding='utf-8')
    for K in result:
        t.write(str(K) + '\n')
    t.close()

def decitoBin(numb):
    # checking if the given number is greater than 1
    if numb > 1:
        # if it is greater than 1 then use recursive approach by dividing number by 2
        decitoBin(numb // 2)
    # printing the binary representation of the given number
    print(numb % 2, end='')

results_path = 'raw_results/GPU_11_15.out'
counter = 0
results = read_file(results_path)
N=len(results)
start_global = timeit.default_timer()
start = timeit.default_timer()
list_of_seeds = []
for k in range(N):
    K = json.loads(results[k])
    K_sp = sc.PureSimplicialComplex(K)
    if K_sp.Pic==4 and K_sp.is_a_seed():
        counter+=1
        list_of_seeds.append(K)
    if k%1000 == 0 and k>1:
        print("Pourcentage effectué: ",(k/N)*100,'%'," Fréquence de seed: ", (counter/k)/100,'% ',"Temps passé: ",timeit.default_timer()-start)
        start = timeit.default_timer()
text(list_of_seeds)
print("Temps passé: ",timeit.default_timer()-start_global)