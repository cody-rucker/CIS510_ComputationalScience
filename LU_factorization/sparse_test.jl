using LinearAlgebra
using SparseArrays

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
F = lu(A)

t1 = @elapsed y = F.L \ b

t2 = @elapsed x = F.U \ y

print(t1)
print("\n")
print(t2)
print("\n")

A2 = spzeros(N,N)

for i = 2:N-1
    A2[i,i] = -2
    A2[i, i-1] = 1
    A2[i, i+1] = 1
end
A2[1, 1] = -2
A2[1, 2] = 1
A2[N, N] = -2
A2[N, N-1] = 1

# LU decomposition using native julia function
F = lu(A2)

t3 = @elapsed y = F.L \ b

t4 = @elapsed x = F.U \ y

print(t3)
print("\n")
print(t4)
print("\n")

print(sizeof(A))
print("\n")
print(sizeof(A2))
