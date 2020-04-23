using LinearAlgebra
using SparseArrays

function comp_piv_matrices(p,q,r)

    N = length(p)

    Id = sparse(1:N,1:N,ones(N),N,N)
    P = copy(Id)
    Q = copy(Id)

    for i = 1:N
        P[i,:] = Id[p[i],:]
        Q[:,i] = Id[:,q[i]]
    end

    R = sparse(1:N,1:N,r,N,N)

    return (P,Q,R)

end

N = 100
A = zeros(N, N)
b = ones(N, 1)

for i = 2:N-1
    A[i,i] = -2
    A[i, i-1] = 1
    A[i, i+1] = 1
end
A[1, 1] = -2
A[1, 2] = 1
A[N, N] = -2
A[N, N-1] = 1

# LU decomposition using native julia function
@time F = lu(A)

@time x = F.U \ (F.L \ b)
@assert A*x ≈ b

# Sparse A

Id = 1:N
Jd = copy(Id)
Vd = -2 * ones(N)

Adiag = sparse(Id, Jd, Vd, N, N)

Il = 2:N
Jl = 1:N-1
Vl = ones(N-1)

Iu = 1:N-1
Ju = 2:N
Vu = ones(N-1)

Asparse = Adiag + sparse(Il, Jl, Vl, N, N) + sparse(Iu, Ju, Vu, N, N)

# LU decomposition using native julia function
@time F = lu(Asparse)

(P,Q,R) = comp_piv_matrices(F.p,F.q,F.Rs)

#time actual solve
@time x = Q*(F.U\(F.L\(P*R*b)))

#time factorization plus solve together:
xMod = Asparse\b

@assert xMod ≈ x

@show sizeof(Asparse)
@show sizeof(A)
