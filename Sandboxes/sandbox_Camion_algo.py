import SimplicialComplex as sc
import numpy as np

A = np.array([[1,1,0,0],
              [0,1,1,0],
              [0,0,1,1],
              [1,0,0,1],
              [1,1,1,0],
              [1,1,0,1],
              [1,0,1,1],
              [0,1,1,1],
              [1,1,1,1],
              ])

B = np.array([[1,1,0,0],
              [0,1,1,0],
              [0,0,1,1],
              [1,0,0,1],
              ])

list_cycles = []


K = sc.PureSimplicialComplex([[1,2,3],[1,2,6],[1,5,3],[1,5,6],[2,3,4],[2,4,6],[3,4,7],[3,5,7],[4,5,7],[4,5,6]])
char_map = sc.Garrison_Scott(K)[0]
relabelling = sc.lifting_algo(K.facets,sc.char_funct_bin_to_numpy(char_map,K.n))
print(relabelling)