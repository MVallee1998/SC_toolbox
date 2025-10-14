using Oscar

function all_nonzero_binary_vectors(n)
    k = 2^n - 1
    mat = Matrix{Int}(undef, k, n)
    row = 1
    for i in 1:k
        v = digits(i, base=2, pad=n)
        mat[row, :] = reverse(v)
        row += 1
    end
    return mat
end

global A = matrix(GF(2), all_nonzero_binary_vectors(5))
global M0 = matroid_from_matrix_rows(A)

# à corriger car problème de rang...

function find_lower_dim_matroids(list_bin_mat)
    local S_rep = Set{Matroid}()
    for M in list_bin_mat
        for v in matroid_groundset(M)
            if [v] in cocircuits(M) # si v est un loop, on skip
                continue
            end
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

global m=31
global S = Set{Matroid}()
push!(S,M0)
while m > 5
    io=open("bin_mat_" * string(m) * "_" * string(m-4),"w")
    print("number of elements ",m," number of matroids ",length(S),"\n")
    for M in S
        println(io,cobases(M))
    end
    close(io)
    global S = find_lower_dim_matroids(S)
    global m-=1
end




