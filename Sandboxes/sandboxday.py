import SimplicialComplex as sc

K=sc.PureSimplicialComplex(None,[[1,5],[2,6],[3,7],[4,8]],4)
K.compute_facets_from_MNF_set()
print(K.facets_bin)