using LinearAlgebra
using Printf

"""
    find_pivot(a)

Locate the element of maximal absolute value in the one-dimensional array a.

# Examples
```jldoctest
julia> A = [ 1 2 -3]
julia> j = find_pivot(A)
```
"""
function find_pivot(a, k)
    N = length(a)
    p = 0

    a = broadcast(abs, a)

    for i = k:N
        if abs(p) < abs(a[i])
            p = a[i]
        end
    end

    q = findall(x -> x==p, a)
    q = q[1]

    return q
end

"""
    computeLUP(a)

Compute and return LUP factorization of square matrix a.

# Examples
```jldoctest
julia> A = rand(3,3)
julia> (L, U, P) = LUPsolve(A)
```
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


    j = find_pivot(Atilde[:,1], 1)

    if j != 1
        (P_ij[j,: ], P_ij[1, :]) = ( P_ij[1, :], P_ij[j,: ] )
        (Atilde[j, 1:N ],Atilde[1, 1:N]) = ( Atilde[1, 1:N], Atilde[j, 1:N ] )
    end

    for k = 1:N-1  # marching across columns


        ell .= Id
        ell_inv .= Id


        for i = k+1:N
            ell[i,k] = -Atilde[i,k] / Atilde[k,k]
            ell_inv[i,k] = Atilde[i,k] / Atilde[k,k]
        end

        Atilde .= ell * Atilde
        L .= L * ell_inv

        j = find_pivot(Atilde[:,k+1], k+1)

        if j != k+1
            (P_ij[j,: ], P_ij[k+1, :]) = ( P_ij[k+1, :], P_ij[j,: ] )
            (L[k+1,1:k ], L[j, 1:k]) = ( L[j, 1:k], L[k+1,1:k ] )
            (Atilde[j, : ],Atilde[k+1, :]) = ( Atilde[k+1, :], Atilde[j, : ] )
        end



    end
    U = Atilde
    P = P_ij

    return (L, U, P)
end

"""
    LUPsolve(A, b)

Solve the matrix equation Ax=b using LU factorization.

# Examples
```jldoctest
julia> A = rand(3,3)
julia> b = rand(3,1)
julia> x = LUPsolve(A, b)
```
"""
function LUPsolve(A, b)

    N = size(L)[1]

    L, U, P = computeLUP(A)

    b = P*b

    y = forward_sub(L, b)
    x = backward_sub(U, y)

    return x



end

"""
    backward_sub(A, b)

Solve the matrix equation Ax=b for upper triangluar matrix A using back substitution.

# Examples
```jldoctest
julia> A = [1 2 3;
            0 4 5;
            0 0 6]
julia> x = backward_sub(A, b)
```
"""
function backward_sub(A, b)
    N = size(A)[2]
    for i = 1:N # march along columns
        k = N+1-i # backwards

        b[k] = b[k] / A[k,k] # normalize the pivot
        for j = 1:k-1

            b[j] = b[j] - b[k] * A[j,k]
        end
    end
    return b
end

"""
    forward_sub(A, b)

Solve the matrix equation Ax=b for lower triangluar matrix A using forward substitution.

# Examples
```jldoctest
julia> A = [1 0 0;
            2 3 0;
            4 5 6]
julia> x = forward_sub(A, b)
```
"""
function forward_sub(A, b)
    N = size(A)[2]

    for i = 1:N # march along columns

        b[i] = b[i] / A[i,i]
        for j = i+1:N

            b[j] = b[j] - b[i]*A[j,i]
        end
    end
    return b
end


#=
N= 10
A = rand(N,N)
b = rand(N,1)

x = LUPsolve(A,b)
x2 = A\b

norm(x-x2)


# displaying run times
N = 10
A = Array{Float64}(undef,N,N)
A .= rand(N,N)

computeLUP(A)
@time (myL, myU, myP) = computeLUP(A)
@assert myL*myU ≈myP*A

b = rand(N,1)
LUPsolve(myL, myU, myP, b)
@time x = LUPsolve(myL, myU, myP, b)
@assert A*x ≈ b
=#



"""
    conj_grad(A, x0, b, tol, iter_max))

Solve the matrix equation Ax=b using iterated conjugate gradient method.

# Arguments
- `tol::Real`: the desired tolerance on relative error.
- `iter_max::Integer`: the maximum number of allowed iterations.

# Examples
```jldoctest
julia> A = rand(3,3)
julia> b = rand(3,1)
julia> x0 = [0; 0; 0]
julia> x = conj_grad(A, x0, b, tol, iter_max))
```
"""
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

    @assert (norm( A*x -b) / norm(x)) <= tol

    return x
end
