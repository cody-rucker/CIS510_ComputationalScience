using LinearAlgebra
using Printf

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

    for k = 1:N-1  # marching across columns

        ell .= Id
        ell_inv .= Id
        for i = k+1:N
            ell[i,k] = -Atilde[i,k] / Atilde[k,k] # compute elimination factors
            ell_inv[i,k] = Atilde[i,k] / Atilde[k,k]
        end


        Atilde .= ell * Atilde
        L      .= L * ell_inv



    end
    U = Atilde

    return (L, U)
end



N = 300
A = Array{Float64}(undef,N,N)
A .= rand(N,N)#[6 -2 2;12 -8 6;3 -13 3]

(myL, myU) = computeLUP(A)
@assert myL*myU ≈ A

b = rand(N,1) # defines the right hand side of Ax = b

#x = LUPsolve(myL, myU, b)

println("Compute my LU factorization")
@time (myL, myU) = computeLUP(A)
@printf "norm(A-LU) = \x1b[31m %e \x1b[0m\n" norm(myL*myU-A)
println("-----------")
println()

#=
println("Peform my LU solve")
@time x = LUPsolve(myL, myU, b)
@printf "norm(Ax-b) = \x1b[31m %e \x1b[0m\n" norm(A*x-b)
println("-----------")
println()
=#



function find_pivot(a)
    k = length(a)
    p = 0

    for i = 1:k
        if p < a[i]
            p = a[i]
        end
    end

    q = findall(x -> x==p, a)

    return p, q
end
