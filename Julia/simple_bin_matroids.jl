using Oscar
using IterTools
using ProgressBars
using LinearAlgebra

function all_nonzero_binary_vectors(n)
    k = (2^n) - 1
    mat = Matrix{Int}(undef, k, n)
    for i in 1:k
        v = digits(i, base=2, pad=n)
        mat[i,:]= v
    end
    return mat
end



# à corriger car problème de rang...

function find_lower_dim_matroids(list_bin_mat)
    local S_rep = Set{Matroid}()
    for M in list_bin_mat
        for v in matroid_groundset(M)
            if v in coloops(M)
                continue
            end
            M1 = deletion(M,v)
            if rank(M)!=rank(M1)
                continue
            end
            is_isom=false
            for M2 in S_rep
                if is_isomorphic(M1,M2)
                    is_isom = true
                    break
                end
            end
            if is_isom==false
                
                push!(S_rep,M1)
            end
        end
    end
    return S_rep
end

global m=31
global A = matrix(GF(2), all_nonzero_binary_vectors(5))
global M0 = matroid_from_matrix_rows(A)
global S = Set{Matroid}()
push!(S,M0)
while m > 5
    io=open("bin_mat_better_" * string(m) * "_" * string(m-5),"w")
    print("number of elements ",m," number of matroids ",length(S),"\n")
    for M in S
        println(io,cobases(M))
    end
    close(io)
    global S = find_lower_dim_matroids(S)
    global m-=1
end

# for i=2:m
#     S=Set{Matroid}()
#     for J in ProgressBar(IterTools.subsets(axes(A,1),i+5))
#         M1=matroid_from_matrix_rows(A[J,:])
#         is_isom=false
#         for M2 in S
#             if is_isomorphic(M1,M2)
#                 is_isom=true
#                 break
#             end
#         end
#         if !is_isom
#             push!(S,M1)
#         end
#     end
#     io=open("bin_mat_" * string(i+5) * "_" * string(i),"w")
#     print("number of elements ",m," number of matroids ",length(S),"\n")
#     for M in S
#         println(io,cobases(M))
#     end
#     close(io)
# end






