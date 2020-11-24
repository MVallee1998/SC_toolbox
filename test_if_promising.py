import Partial_Enum_PLS_IDCM
import SimplicialComplex as sc
import json

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

results = read_file('result/PLS_test')

K_result=[]
l=0

for K_bytes in results:
    K_bin= json.loads(K_bytes)
    K_sp = sc.PureSimplicialComplex(K_bin)
    if K_sp.is_Z2_homology_sphere():
        if K_bin not in K_result:
            K_result.append(K_bin)

    # list_2_pow = [1]
    # for k in range(8):
    #     list_2_pow.append(list_2_pow[-1]*2)
    # K= []
    # for facet_bin in K_bin:
    #     K.append([])
    #     for i in range(len(list_2_pow)):
    #         element = list_2_pow[i]
    #         if element | facet_bin == facet_bin:
    #             K[-1].append(i+1)
    # K.sort()
    # Chain, Ind = Partial_Enum_PLS_IDCM.Chain_cpx(K, 6)
    # if Partial_Enum_PLS_IDCM.is_Homology_Sphere(Chain, Ind) and K not in K_result:
    #     K_result.append(K)
print(K_result)