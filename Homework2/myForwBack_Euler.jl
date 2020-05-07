using LinearAlgebra
using Plots

function my_forward_Euler!(Y, Δt, t, A, F, y, x, N, M)
    #M = length(t)
    #t = t1:Δt:tf
    #M = Integer(ceil((tf-t1)/Δt))
    Y[2:N,1] = y[:]

    for n = 1:M
        b = Δt * F(x[2:N],t[n])
        y[:] = A * y[:] + b
        Y[2:N,n] .= y[:]

    end

    #return (t, Y)

end


function my_backward_Euler!(Y, Δt, t, A, F, y, x, N, M)


    Y[2:N,1] = y[:]

    #Exact = Matrix{Float64}(undef,2,N+1)
    #Exact[:,1] = y[:]

    Id = I(N-1) #Matrix{Float64}(I,N+1,N+1)

    for n = 2:M-1
        y[:] = A\(y[:] .+ Δt*F(x[2:N], t[n]))
        Y[2:N,n] = y[:]
        #Exact[:,n] = exact(t[n])
    end
end

#=
A = [-3.0 13.0;-5.0 -1.0]
y1 = Array{Float64}(undef,2,1)
y1 = [3.0;-10.0]
Δt = 0.01
t1 = 0
tf = 3
t = t1:Δt:tf
N = Integer((tf - t1)/Δt)  # total number of points is N + 1
c = zeros(length(t))

Y = Matrix{Float64}(undef,2,N+1)
#(t, Exact) = my_backward_Euler!(Y, Δt, t1, tf, A, c, y1)
(t, Exact) = my_forward_Euler!(Y, Δt, t1, tf, A, c, y1)

labels = ["approx_y1" "exact_y1" "approx_y2" "exact_y2"]
markershapes = [:circle :star :circle :star]
markercolors = [:green :green :purple :purple]
strd = 1


#plot(t[1:strd:end],[Y[1,1:strd:end], Exact[1,1:strd:end], Y[2,1:strd:end], Exact[2,1:strd:end]], label = labels, shape = markershapes, color = markercolors)
plot(Y[1,1:strd:end],Y[2,1:strd:end])
=#
