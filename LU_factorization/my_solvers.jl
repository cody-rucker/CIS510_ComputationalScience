using LinearAlgebra
using Printf


function find_pivot(a)
    k = length(a)
    p = 0

    a = broadcast(abs, a)

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


        ell = copy(Id)

        #j = find_pivot(A[:, k])
        j = findmax(Atilde[k:N ,k])[2]


        if j != k
            (P_ij[j,: ], P_ij[k, :]) = ( P_ij[k, :], P_ij[j,: ] )
            #(ell[k,1:k-1 ], ell[j, 1:k-1]) = ( ell[j, 1:k-1], ell[k,1:k-1 ] )
            (Atilde[j, k:N ],Atilde[k, k:N]) = ( Atilde[k, k:N], Atilde[j, k:N ] )
        end

        @assert Atilde[k,k] != 0


        #P .= P_ij * P

        for i = k+1:N
            ell[i,k] = -Atilde[i,k] / Atilde[k,k] # compute elimination factors
            Atilde[i, k:N] = Atilde[i, k:N]  + ell[i,k] * Atilde[k, k:N]
            ell_inv[i,k] = Atilde[i,k] / Atilde[k,k]
        end


        #Atilde .= ell * Atilde
        #L      .=  ell_inv * L
        L[:, k] = ell[:, k]
        #L .= L * ell_inv


    end
    U = Atilde
    P = P_ij

    return (L, U, P)
end



N = 5
A = Array{Float64}(undef,N,N)
A .= rand(N,N)#[6 -2 2;12 -8 6;3 -13 3]

(myL, myU, myP) = computeLUP(A)
#@assert myL*myU*myP ≈ A

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

#=

function conj_grad(A, x0, b, tol, iter_max)

    i = 0
    x = x0
    r = b - (A * x)
    d = r
    δ_new = r' * r
    δ0 = δ_new

    while i < (iter_max) && δ_new > (tol)^2 * δ0

        q = A * d
        α = δ_new / (d' * q)
        x = x + α * d

        # r = r-αq is less computationally expensive than r = b - Ax
        # but r = r - αq is generated without any feedback from the
        # value of x_i. So the accumulation of floating point roundoff
        # error may cause x_i to converge to some point near x.
        # By periodically using r = b - Ax to recompute the correct
        # residual we can avoid this issue of accumulated roundoff error.
        if  i%50 == 0
            r = b - A * x
        else
            r = r - α * q
        end

        δ_old = δ_new
        δ_new = r' * r
        β = δ_new / δ_old
        d = r + β * d

        i = i+1


        # if err_R <= tol we are within our desired tolerence for error
        # and can terminate the algorithm early.
        if (norm( A*x -b) / norm(x)) <= tol
            break
        end
    end


    return x
end

N = 10
B = rand(N, N)
A = B' * B

x0 = vec(zeros(N,1))
b = vec(ones(N,1))
tol = 10^-6
iter_max = 100

x = conj_grad(A, x0, b, tol, iter_max)
@assert (norm( A*x -b) / norm(x)) <= tol
=#
