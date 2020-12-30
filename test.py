import SimplicialComplex as sc
import numpy as np
import timeit
import graph_methods as gm
import linear_alg_method as lam

# char_funct = [12, 7, 11, 13, 14, 15, 8, 4, 2, 1]
# m = 10
# n = 6
# M, facets, ridges = lam.construct_matrix(char_funct, n, m)
# K = []
#
# for face in facets:
#     K.append(sc.binary_to_face(face, m))

K_sp = sc.PureSimplicialComplex([31, 55, 79, 91, 94, 103, 115, 118, 157, 181, 205, 217, 220, 229, 241, 244])
K_sp.create_FP()
print(K_sp.is_Z2_homology_sphere(),K_sp.is_promising(),K_sp.is_closed())