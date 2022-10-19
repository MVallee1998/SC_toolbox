import SimplicialComplex as sc
import Z2_linear_algebra as Z2la
import Puzzle_algorithm as pa
import Z2_Cohomology as Z2c
from sympy import *
import numpy as np
from sympy import groebner
x, y, z = symbols('x y z')
from sympy.abc import x, y, z

import json
P_5 = sc.PureSimplicialComplex([[1,2],[1,5],[2,3],[3,4],[4,5]])

P = sc.multiple_wedge(P_5,[1,1,1,0,0])
Q = sc.multiple_wedge(P_5,[1,1,1,0,0])

GL3 = []
for i0 in range(1, 8):
    for i1 in range(1, 8):
        for i2 in range(1, 8):
                A_bin = Z2la.Z2Array(3, [i0, i1, i2])
                if A_bin.is_invertible():
                    A = np.zeros((3, 3))
                    A[0] = sc.int_to_bin_array(i0, 3)
                    A[1] = sc.int_to_bin_array(i1, 3)
                    A[2] = sc.int_to_bin_array(i2, 3)
                    GL3.append(A.copy())

base = [Poly(x, gens=Z2c.list_gens, domain=GF(2)),Poly(y, gens=Z2c.list_gens, domain=GF(2)),Poly(z, gens=Z2c.list_gens, domain=GF(2))]
def transform_base(base,G):
    new_base = []
    for i in range(len(base)):
        for j in range((len(base))):
            if G[i,j]==1:
                if len(new_base)<i+1:
                    new_base.append(base[j])
                else:
                    new_base[-1]+=base[j]
    return new_base


list_IDCM_P = [sc.DCM_bin_to_IDCM_bin(DCM_bin,P.n) for DCM_bin in sc.Garrison_Scott(P)]
# list_IDCM_Q = [sc.DCM_bin_to_IDCM_bin(DCM_bin,Q.n) for DCM_bin in sc.Garrison_Scott(Q)]

# list_IDCM_P = sc.IDCM_Garrison_Scott(P_5)

isom_cohom_P = []
for i in range(len(list_IDCM_P)):
    print("global: ",(i/len(list_IDCM_P))*100,'%')
    isom=False
    IDCM_2 = list_IDCM_P[i]
    for IDCM_1 in isom_cohom_P:
        GB_1 = Z2c.construct_ideal(P,IDCM_1,base)
        for G in GL3:
            GB_2 = Z2c.construct_ideal(P,IDCM_2,transform_base(base,G))
            if GB_1==GB_2:
                isom=True
                break
        if isom:
            break
    if not isom:
        isom_cohom_P.append(IDCM_2.copy())
print(len(isom_cohom_P))

# isom_cohom_Q = []
# for i in range(len(list_IDCM_Q)):
#     print("global: ",(i/len(list_IDCM_Q))*100,'%')
#     isom=False
#     IDCM_2 = list_IDCM_Q[i]
#     for IDCM_1 in isom_cohom_Q:
#         GB_1 = Z2c.construct_ideal(Q,IDCM_1,base)
#         for G in GL3:
#             GB_2 = Z2c.construct_ideal(Q,IDCM_2,transform_base(base,G))
#             if GB_1==GB_2:
#                 isom=True
#                 break
#         if isom:
#             break
#     if not isom:
#         isom_cohom_Q.append(IDCM_2.copy())
# print(len(isom_cohom_Q))
#
# for i in range(len(isom_cohom_P)):
#     print("global: ",(i/len(isom_cohom_P))*100,'%')
#     IDCM_P = isom_cohom_P[i]
#     GB_P = Z2c.construct_ideal(P,IDCM_P,base)
#     for j in range(len(isom_cohom_Q)):
#         IDCM_Q=isom_cohom_Q[j]
#         for G in GL3:
#             GB_Q = Z2c.construct_ideal(Q,IDCM_Q,transform_base(base,G))
#             if GB_Q==GB_P:
#                 print('hello')
#                 break