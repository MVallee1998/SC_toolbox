for a in range(10):
    for b in range(10):
        for c in range(10):
            for d in range(10):
                P = Polyhedron(ieqs=[(0,1,0,0,0,0,0),(0,0,1,0,0,0,0),(0,0,0,1,0,0,0),(0,0,0,0,1,0,0),(0,0,0,0,0,1,0),(a,-1,-1,0,-1,0,0),(b,-1,0,-1,0,-1,0),(c,0,-1,0,-1,-1,0),(d,0,0,-1,-1,-1,0)],eqns=[(-1,0,0,0,0,0,1)])
                if P.is_lattice_polytope() and P.is_simple() and len(P.f_vector())==7 and P.f_vector()[5]==9:
                    print(P.f_vector())
                    C = Cone(P.Vrepresentation())
                    for v in C.Hilbert_basis():
                        if v[5]!=1:
                            print('hello')