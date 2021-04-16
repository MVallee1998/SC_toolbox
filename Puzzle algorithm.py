import SimplicialComplex as sc
import numpy as np
import json
import linear_alg_method as lam
import timeit

list_2_pow = [1]
for k in range(16):
    list_2_pow.append(list_2_pow[-1] * 2)


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data




def puzzle_algo(K,J):
    p = K.Pic
    n = K.n
    m = K.m
    list_IDCM_bin = sc.IDCM_Garrison_Scott(K)
    list_IDCM= []
    for IDCM_bin in list_IDCM_bin:
        IDCM = np.zeros((m,p))
        for i in range(m):
            for j in range(p):
                if list_2_pow[j] | IDCM_bin[i] == IDCM_bin[i]:
                    IDCM[i,j] = 1
        list_IDCM.append(IDCM.copy())
    list_CM_bin = sc.Garrison_Scott(K)
    list_CM_1 = []
    for CM_bin in list_CM_bin:
        CM = np.zeros((n,m))
        for j in range(m):
            for i in range(n):
                if CM_bin[j] | list_2_pow[i] == CM_bin[j]:
                    CM[i,j] = 1
        list_CM_1.append(CM.copy())
    list_CM= []
    for IDCM in list_IDCM:
        CM = np.zeros((n,m))
        CM[:,:n] = np.eye(n)
        CM[:,n:m] = IDCM[:n,:]
        list_CM.append(CM.copy())
    for CM1 in list_CM_1:
        is_in = False
        for CM in list_CM:
            if (CM1==CM).all():
                is_in = True
                break
        if not is_in:
            print("coucou")


for n in range(2,8):
    m=n+4
    results = read_file('final_results_BAK/PLS_%d_%d' % (m, n))
    start = timeit.default_timer()
    for K_byte in results:
        K = sc.PureSimplicialComplex(json.loads(K_byte))
        puzzle_algo(K,1)



