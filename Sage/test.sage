GL4 = GL(4,GF(2))
P.<x,y,z,t> = PolynomialRing (GF(2) ,order ='degrevlex')
import numpy as np

list_poly_gens = [x,y,z,t]


print(len(GL4))
def transform_base(base, G):
    new_base = []
    for i in range(len(base)):
        for j in range((len(base))):
            if G[i][j] == 1:
                if len(new_base) < i + 1:
                    new_base.append(base[j])
                else:
                    new_base[-1] += base[j]
    return new_base


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


def construct_ideal(MNF_set_bin, IDCM_bin, test,m,n):
    p = m - n
    list_monomials = []
    for k in range(n):
        P = 0
        for l in range(p):
            if (IDCM_bin[k] >> l) & 1:
                P += test[p-l-1]
        list_monomials.append(P)
    list_monomials += test
    # print(list_monomials)
    L = []
    for MNF_bin in MNF_set_bin:
        P = 1
        for l in range(m):
            if (MNF_bin >> l) & 1:
                P *= list_monomials[l]
        L.append(P)
    L_ideal = Ideal(L)
    GB = L_ideal.groebner_basis()
    return GB

Hexagon_MNF = [5, 9, 10, 17, 18, 20, 34, 36, 40]
IDCM_1=[10,5,8,4,2,1]
IDCM_2=[11,5,8,4,2,1]
GB_1 = construct_ideal(Hexagon_MNF,IDCM_1,list_poly_gens,6,2)
for k in range(len(GL4)):
    if k%1000==0:
        print((k/len(GL4))*100,'%')
    G=GL4[k].list()
    new_base=transform_base(list_poly_gens, G)
    GB_2 = construct_ideal(Hexagon_MNF, IDCM_2, new_base,6,2)
    if GB_1 == GB_2:
        print(G)
        isom = True
        break
print("finished")