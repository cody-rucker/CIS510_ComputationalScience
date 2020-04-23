using LinearAlgebra
using IterativeSolvers

N = 1000
B = rand(N,N)
A = B*B'

b = ones(N,1)

# LU solve
F = lu(A)

@time x = F.U \ (F.L \ b)
@assert A*x ≈ b

# CG solve
@time y = cg(A, b; tol=1e-9, maxiter=3*N)
@assert A*y ≈ b
