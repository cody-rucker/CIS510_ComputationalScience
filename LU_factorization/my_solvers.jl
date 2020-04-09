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
      ℓ   = copy(Id)
      ℓ⁻¹ = copy(Id)
      Ã   = copy(A)
      L   = copy(Id)


      for k = 1:N-1
            ℓ   .= Id
            ℓ⁻¹ .= Id


            for i = k+1:N
                  ℓ[i, k]   = -Ã[i,k] / Ã[k,k]
                  ℓ⁻¹[i, k] = Ã[i,k] / Ã[k,k]
            end

            Ã .= ℓ * Ã
            L .= L * ℓ⁻¹
      end

      U = Ã

      return L, U
end


N = 100
A = Array{Float64}(undef, N, N)
A .= rand(N, N)

(myL, myU) = ComputeLU(A)

@assert myL*myU ≈ A
