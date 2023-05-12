import SimplicialComplex as sc
import Z2_linear_algebra as Z2la
import Puzzle_algorithm as pa
import Z2_Cohomology as Z2c
from sympy import *
import numpy as np

Hexagon = sc.PureSimplicialComplex([[1,2],[1,6],[2,3],[3,4],[4,5],[5,6]])
Hexagon.compute_MNF_set()
print(Hexagon.MNF_set_bin)
IDCM_1=[10,5,8,4,2,1]
IDCM_2=[11,5,8,4,2,1]

GB_1 = Z2c.construct_ideal(Hexagon, IDCM_1, Z2c.list_poly_gens)
# GB_2 = Z2c.construct_ideal(Hexagon, IDCM_2, Z2c.list_poly_gens)

# for k in range(len(Z2c.GL4)):
#     if k%1000==0:
#         print((k/len(Z2c.GL4))*100,'%')
#     G=Z2c.GL4[k]
#     new_base=Z2c.transform_base(Z2c.list_poly_gens, G)
#     GB_2 = Z2c.construct_ideal(Hexagon, IDCM_2, new_base)
#     if GB_1 == GB_2:
#         print(G)
#         isom = True
#         break
# print("finished")
