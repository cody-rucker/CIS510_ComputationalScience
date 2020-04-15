using LinearAlgebra
using Printf


function find_pivot(a)
    k = length(a)
    p = 0

    for i = 1:k
        if abs(p) < abs(a[i])
            p = a[i]
        end
    end

    q = findall(x -> x==p, a)
    q = q[1]

    return q
end

"""
    LUPsolve(a)

Compute and return LU factorization of square matrix a.

# Examples
'''
julia> A = rand(3,3)
julia> (L, U) = LUPsolve(A)

etc

'''
"""
function computeLUP(A)
    # TODO: develop to include partial pivoting
    N = size(A)[1]

    Id      = Matrix{Float64}(I, N, N)
    ell     = copy(Id)
    ell_inv = copy(Id)
    Atilde  = copy(A)
    L       = copy(Id)
    P       = copy(Id)
    P_ij    = copy(Id)

    for k = 1:N-1  # marching across columns

        ell .= Id
        ell_inv .= Id

        j = find_pivot(A[:, k])
        #j = findmax(Atilde[:,k])[2]


        if j != k
            (P_ij[j,: ], P_ij[k, :]) = ( P_ij[k, :], P_ij[j,: ] )
            (ell[k,1:k-1 ], ell[j, 1:k-1]) = ( ell[j, 1:k-1], ell[k,1:k-1 ] )
            (Atilde[j,: ],Atilde[k, :]) = ( Atilde[k, :], Atilde[j,: ] )
        end


        #Atilde .=P_ij * Atilde
        print(Atilde[k,k])
        print("\n")
        @assert Atilde[k,k] != 0


        #P .= P_ij * P

        for i = k+1:N
            ell[i,k] = Atilde[i,k] / Atilde[k,k] # compute elimination factors
            Atilde[i, k:N] = Atilde[i, k:N] - ell[i,k] * Atilde[k, k:N]
            ell_inv[i,k] = Atilde[i,k] / Atilde[k,k]
        end


        #Atilde .= ell * Atilde
        #L      .=  ell_inv * L
        L[:, k] = ell[:, k]




    end
    U = Atilde
    P = P_ij

    return (L, U, P)
end



N = 5
A = Array{Float64}(undef,N,N)
A .= rand(N,N)#[6 -2 2;12 -8 6;3 -13 3]

(myL, myU, myP) = computeLUP(A)
@assert myL*myU*myP ≈ A

b = rand(N,1) # defines the right hand side of Ax = b

#x = LUPsolve(myL, myU, b)

println("Compute my LU factorization")
@time (myL, myU, myP) = computeLUP(A)
@printf "norm(A-LU) = \x1b[31m %e \x1b[0m\n" norm(myL*myU-myP*A)
println("-----------")
println()

#=
println("Peform my LU solve")
@time x = LUPsolve(myL, myU, b)
@printf "norm(Ax-b) = \x1b[31m %e \x1b[0m\n" norm(A*x-b)
println("-----------")
println()
=#
