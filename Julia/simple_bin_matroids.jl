using Oscar

global A = matrix(GF(2),[[1,1,0,0],[1,0,1,0],[1,0,0,1],[0,1,1,0],[0,1,0,1],[0,0,1,1],[1,1,1,0],[1,1,0,1],[1,0,1,1],[0,1,1,1],[1,1,1,1],[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
global M0 = matroid_from_matrix_rows(A)

# à corriger car problème de rang...

function find_lower_dim_matroids(list_bin_mat)
    local S_rep = Set{Matroid}()
    for M in list_bin_mat
        for v in matroid_groundset(M)
            M1 = deletion(M,v)
            is_isom=false
            for M2 in S_rep
                if is_isomorphic(M1,M2)
                    is_isom = true
                    break
                end
            end
            if is_isom==false
                str = revlex_basis_encoding(M1)
                push!(S_rep,matroid_from_revlex_basis_encoding(str,rank(M1),length(M1)))
            end
        end
    end
    return S_rep
end

global m=15
global S = Set{Matroid}()
push!(S,M0)
while m > 5
    io=open("bin_mat_" * string(m) * "_" * string(m-4),"w")
    for M in S
        println(io,cobases(M))
    end
    close(io)
    global S = find_lower_dim_matroids(S)
    global m-=1
end




