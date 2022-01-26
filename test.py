import SimplicialComplex as sc
import json
import timeit
import numpy as np
import collections

m = 15
n = 11


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def text(result):
    name = 'result/PLS_%d_%d_good_seeds_lin_alg' % (m, n)
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

results = read_file('partial_results/PLS_%d_%d' % (m, n))
l=0
counter = 0
start = timeit.default_timer()
for K_bytes in results:
    l += 1
    K = json.loads(K_bytes)
    K_sp = sc.PureSimplicialComplex(K)
    K_sp.compute_MNF_set()
    K_sp.MNF_bin_to_MNF()
    counter= 0
    facets_test = []
    for v in K_sp.facets_bin:
        facets_test.append(sc.binary_to_face(v ^ (sc.list_2_pow[15]-1),15))
    K_sp.create_f_vector()
    for IDCM in sc.IDCM_Garrison_Scott(K_sp)[1:]:
        char_map = np.zeros((m,m-n),dtype=int)
        for k in range(m):
            for i in sc.binary_to_face(IDCM[k],m):
                char_map[k,i-1] = 1
        print(char_map)
        # char_map[8,1]=2
        # char_map[3, 1] = 2
        char_map[4,1]=-1
        char_map[10,2]=2
        print(char_map)
        facets_test_np = np.array(facets_test,dtype=int)
        facets_test_np-=1
        for facet in facets_test_np:
            if np.abs(np.linalg.det(char_map[facet,:])) !=1:
                print("this one",facet+1,'\n',char_map[facet,:],np.abs(np.linalg.det(char_map[facet,:])))
    print(counter,K_sp.f_vector,K_sp.facets)
    print(K_sp.MNF_set)
