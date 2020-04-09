using LinearAlgebra

"""
    ComputeLU(a)

Compute and return LU factorization of sqaure matrix a.

# Examples
'''
julia> A = rand(3, 3)
julia> (L, U) = ComputeLU(A)
'''
"""


function ComputeLU(A)
      N = size(A)[1]

      Id  = Matrix{Float64}(I, N, N)
      ell   = copy(Id)
      ell_inv = copy(Id)
      Atilde   = copy(A)
      L   = copy(Id)


      for k = 1:N-1
            ell   .= Id
            ell_inv .= Id


            for i = k+1:N
                  ell[i, k]   = -Atilde[i,k] / Atilde[k,k]
                  ell_inv[i, k] = Atilde[i,k] / Atilde[k,k]
            end

            Atilde .= ell * Atilde
            L .= L * ell_inv
      end

      U = Atilde

      return L, U
end


N = 100
A = Array{Float64}(undef, N, N)
A .= rand(N, N)

(myL, myU) = ComputeLU(A)

@assert myL*myU ≈ A
